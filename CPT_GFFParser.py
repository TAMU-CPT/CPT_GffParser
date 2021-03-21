# Copyright 2020-2021, Anthony Criscione
# Developed for the Center for Phage Technology, Texas A&M University
#
# Distributed under the BSD 3-Clause License, see included LICENSE file

from Bio import SeqIO, SeqFeature
from Bio.SeqFeature import FeatureLocation, CompoundLocation
from Bio.Seq import Seq, UnknownSeq
from collections import OrderedDict
from Bio.SeqRecord import SeqRecord

from CPT_ParserHelpers import *
from gffSeqFeature import *
import sys

#Try/Except blocks used for limited python 2.7 compatibility. Python3 spec is within the try block
try:
  import urllib.parse
except:
  import urllib

def lineAnalysis(line, codingTypes = ["CDS"]):
    IDName = ""
    startLoc = -1
    endLoc = -1
    scoreIn = 0.0
    if len(line.strip()) == 0 or line.strip() == "\n":
      return None, None, None
    if line[0] == "#":
      if len(line) > 2 and line[1] == "#":
        return None, line[2:-1], None 
      # else handle ## Pragmas
      else: 
        return None, None, None
    errorMessage = ""

    fields = line.split("\t")
    if len(fields) != 9:
      errorMessage += "GFF3 is a 9-column tab-separated format, line has %d columns.\n" % (len(fields))
      if len(fields) > 9:
        errorMessage += "Possible unescaped tab in a qualifier field.\n"
        return errorMessage, None, None

    for x in range(0, len(fields)):
      if fields[x] == "":
        errorMessage += "Field #%d is empty. Please supply correct or default value.\n" % (x + 1)
    if errorMessage != "":
      return errorMessage, None, None

    idEval = validateID(fields[0])
    if len(idEval) != 0:
      errorMessage += "Organism ID contains the following invalid characters: %s\n" % (idEval)

    # fields[1]
    
    # fields[2]

    # fields[3]
    
    if fields[3][0] in "<>":
      uncert = 1
      fields[3] = (fields[3][1:]).strip()
    else:
      uncert = 0
      fields[3] = fields[3].strip()

    if isNum(fields[3]):
      startLoc = int(fields[3])
    else:
      errorMessage += "Feature location start contains non-numeric character.\n"

    # fields[4]
    
    if fields[4][0] in "<>":
      uncert = 1
      fields[4] = (fields[4][1:]).strip()
    else:
      uncert = 0
      fields[4] = fields[4].strip()

    if isNum(fields[4]):
      endLoc = int(fields[4])
    else:
      errorMessage += "Feature location end contains non-numeric character.\n"

    if startLoc >= 0 and endLoc >= 0 and endLoc < startLoc:
      errorMessage += "Feature Location end is less than start (GFF spec requires all features, regardless of strand, to have the lower number as the start).\n"

    # fields[5]
    if fields[5] != ".":
      try:
        scoreIn = float(fields[5])
      except:
        scoreIn = 0.0
        errorMessage += "Score field could not be interpreted as a floating-point (real) number. Ensure notation is correct.\n"

    # fields[6]
    if fields[6] not in ['-', '+', '.', '?']:
      errorMessage += "Feature strand must be '+', '-', '.', or '?', actual value is '%s'.\n" % (fields[6])
    
    # fields[7]
    possibleCodingError = False
    if fields[7] not in ['.', '0', '1', '2']:
      errorMessage += "Expected 0, 1, 2, or . for Phase field value, actual value is '%s'.\n" % (fields[7])
    elif fields[7] =='.' and fields[1] in codingTypes:
      possibleCodingError = True
      # Try to resolve after qualifiers are read in, allow codon_start to supply this value
    if fields[7] == '.':
      phaseIn = 0
    else:
      phaseIn = int(fields[7])

    # fields[8]

    keyName = ""
    valNames = [""]
    valInd = 0
    parseMode = 0  
    qualDict = {} 
    contCounter = 0
    for x in range(0, len(fields[8])):
      currChar = fields[8][x]
      if contCounter:
        contCounter += -1
        continue
      if parseMode == 0:
        if not (currChar in "=,;%\n"):
          keyName += currChar
        elif currChar == "=":
          if len(keyName) == 0:
            errorMessage += "No ID name supplied for a value in the qualifiers field, aborting.\n"
            break
          parseMode = 1
          continue
        elif currChar == "%" and (fields[8][x+1] in encoders) and (fields[8][x+2] in encoders):
          try:
            keyName += urllib.parse.unquote(fields[8][x:x+3])
          except:
            keyName += urllib.unquote(fields[8][x:x+3])
          contCounter = 2
        else: #Encode special char
          encoded = str(hex(ord(currChar)))
          keyName += "%" + encoded[2:].upper()
      elif parseMode == 1:
        if not (currChar in "=,;%\n"):
          valNames[valInd] += currChar
        elif currChar == ",":
          valInd += 1
          valNames.append("")
        elif currChar == "=":
          valNames[valInd] += "%3D"
        elif currChar == "%" and (fields[8][x+1] in encoders) and (fields[8][x+2] in encoders):
          try:
            valNames[valInd] += urllib.parse.unquote(fields[8][x:x+3])
          except:
            valNames[valInd] += urllib.unquote(fields[8][x:x+3])
          contCounter = 2
        elif currChar == "\n":
          parseMode = 2
        else:
          if x == len(fields[8]) - 2: # Assume if last char in fields[8] is a semicolon, then just the end of qualifier 
            parseMode = 2
          elif encodeFromLookahead(fields[8][x+1:]):
            valNames[valInd] += "%3B"
            continue
          else:
            parseMode = 2
          
      if parseMode == 2: # Do not elif, this is used as a wrapup for each qualifier and we want it checked if parsemode == 1 incremented itself 
        if keyName not in qualDict.keys():
          qualDict[keyName] = valNames
        else:
          for x in valNames:
            qualDict[keyName].append(x)
        keyName = ""      
        valNames = [""]
        valInd = 0
        parseMode = 0

    for x in qualDict.keys():
      if x == "ID":
        #if len(qualDict[x]) > 1:
          #errorMessage += "More than one ID supplied for feature.\n"
        IDName = qualDict[x][0]
      if x == "codon_start":
        if fields[7] == ".":
          phaseIn = int(qualDict[x][0])
          possibleCodingError = False
        # Else error? Is overriding a value of 0 desirable?
 
    if possibleCodingError: # No valid value in phase field, no codon_start, and feature is a type listed in codingTypes
      errorMessage += "Expected 0, 1, or 2 in Phase field for %s-type feature, actual value is '%s'.\n" % (fields[1], fields[7])  

    if startLoc == -1 or endLoc == -1 or (not(fields[6] in '-+.?')):
      errorMessage += "Unable to construct feature location, aborting.\n"
    elif fields[6] == '+':
      featLoc = FeatureLocation(startLoc - 1, endLoc, strand = +1)
    elif fields[6] == '-':
      featLoc = FeatureLocation(startLoc - 1, endLoc, strand = -1)
    else:
      featLoc = FeatureLocation(startLoc - 1, endLoc, strand = 0)

    if "Parent" in qualDict.keys():
      for x in qualDict["Parent"]:
        if x == IDName:
          errorMessage += "Feature lists itself as a sub_feature, cycles/loops not permitted in GFF format.\n"
     
    if errorMessage != "":
      return errorMessage, None, None
    return None, fields[0], gffSeqFeature(featLoc, fields[2], '', featLoc.strand, IDName, qualDict, None, None, None, phaseIn, scoreIn, fields[1])   
        
def gffParse(gff3In, base_dict = {}, outStream = sys.stderr, codingTypes=["CDS"], metaTypes = ["remark"], suppressMeta = 2, pragmaPriority = True, pragmaOverridePriority = True):
    # gff3In --- source file
    # base_dict --- file with additional SeqRecord information. Keys are OrganismIDs and values are SeqRecords.
    #               For BCBio backwards compatibility.
    # outStream --- output filestream or stringstream
    # codingTypes --- list of feature types where a non-'.' phase value is expected, passed along to lineAnalysis
    # metaTypes --- list of metadata feature types. Features of this type will be affected by the remaining arguments
    # suppressMeta --- Suppress metadata fields. Int, where 0 == no suppression, all metadata from features and pragmas 
    #                  will be read and output to the SeqRecord as .annotation entries.
    #                  1 == As above, but metadata features will not be entered into the SeqRecord's feature list
    #                  after their metadata is recorded and entered into the SeqRecord.annotation
    #                  2 == Total suppression, no metadata features will even be processed, and no pragmas except
    #                  those related to sequence length (##FASTA and ##sequence-region) and ##gff-version (required 
    #                  by GFF spec) will be utilized.
    # pragmaPriority --- In cases where pragmas and metadata features disagree/conflict, pragmas will take precedence for creating
    #                    SeqRecord.annotation value if true, else the feature will.
    # pragmaOverridePriority --- Similar to above, in the event of a conflict between metadata features and pragmas, the pragma's 
    #                            value will override the metadata gffSeqFeature.qualifier value with the pragma's own. This will force
    #                            the metadata and pragmas to sync, and avoid future discrepancies. Should only be used with pragmaPrority 
 
    if hasattr(gff3In, 'readlines') and callable(getattr(gff3In, 'readlines')):
      gff3In = gff3In.readlines()

    fastaDirective = False # Once true, must assume remainder of file is a FASTA, per spec
    errOut = ""
    warnOut = ""
    lineInd = 0
    pragmaAnnotesDict = {} # Annotations dictionaries, one for ones derived via pragmas and another for ones derived from meta features
    metaAnnotesDict = {}   # Keys are annotation title, values are a list of lists, where the first entry in a list is the value, and 
                           # the rest the orgIDs which correspond to the value
    orgDict = {} # Dictionary of orgIDs, with a list of their features (Top-level or otherwise) as its value
    finalOrg = {}
    seekParentDict = {} # Dictionary of orgIDs, with a list of indexes in the equivalent orgDict list which are subfeatures
    seqDict = {}
    regionDict = {} # Dictionary of max length for records (For handling circular and bound checking)
                    # Values are (maxVal, status), where status is 0 if maxVal is derived from features, 1 if derived from sequence-region 
                    # pragma (Enforce this boundary), or -1 if derived from a feature with is_circular (Unable to construct record)
    currFastaKey = ""
    
    if pragmaPriority:  # In cases where pragmas and meta-features conflict, such as sequence-region
      pragBit = 1       # pragmas will take priority
      annoteBit = 0
    else:               # Else meta features will
      pragBit = 0       # Split into ints to cleanly slot into regionDict format
      annoteBit = 1

    for line in gff3In:
      lineInd += 1
      err = None
      prag = None
      res = None

      ### FASTA Pragma Handling
      if line[0] == ">":	# For compatibility with Artemis-style GFF
        fastaDirective = True
      if not fastaDirective:
        err, prag, res = lineAnalysis(line, codingTypes)
      else:
        if line[0] == ">":
          currFastaKey = line[1:-1]
          if currFastaKey not in seqDict.keys():
            seqDict[currFastaKey] = ""
        elif line[0] == "#":
          continue
        elif line:
          seqDict[currFastaKey] += (line[:-1]).strip()
          continue
      ### Error message construction
      if err:
        errOut += "Line %d: %s\n" % (lineInd, err)
      ### Pragma handling
      if prag and not res:
        prag = prag.split(" ")
        if prag[0] == "FASTA":
          fastaDirective = True
        elif prag[0] == "sequence-region":
          if prag[1] not in regionDict.keys():
            regionDict[prag[1]] = (int(prag[2]) - 1, int(prag[3]), pragBit)
          elif pragBit > regionDict[prag[1]][2]:
            regionDict[prag[1]] = (int(prag[2]) - 1, int(prag[3]), pragBit)
        elif prag[0] == "#":
          orgDict, resolveErr = resolveParent(orgDict, seekParentDict)
          if resolveErr:
            errOut += resolveErr
          finalOrg = rAddDict(finalOrg, orgDict)
          seekParentDict = {}
          orgDict = {}   
        elif suppressMeta < 2: 
          if prag[0] in pragmaAnnotesDict.keys():
            dictVal = " ".join(prag[1:])
            pragmaAnnotesDict[prag[0]].append([dictVal])
          else:
            dictVal = " ".join(prag[1:])
            pragmaAnnotesDict[prag[0]] = [[dictVal]]
      ### Feature Handling
      if res:
        if suppressMeta == 2 and res.type in metaTypes:
          continue
      ## First time encountering orgID
        if prag not in orgDict.keys():
          orgDict[prag] = [res]
          seekParentDict[prag] = []
          possSeq = base_dict.get(prag, None)
          # Process base_dict
          # .seq priority is: ##FASTA directives will always define sequence-region and seq if present (done further down)
          #                   base_dict is next, and will also accept an empty seq, so take care with what's passed in this field
          #                   Finally, parser will infer an UnknownSeq from either ##sequence-region pragma or the 'last' feature,
          #                   depending on arguments passed to parser.
          if isinstance(possSeq, SeqRecord):
            if possSeq.seq == None:
              seqDict[prag] = ""
            else:
              seqDict[prag] = str(possSeq.seq)
          else:
            seqDict[prag] = ""
          for x in pragmaAnnotesDict.keys():
            if prag in pragmaAnnotesDict[x][-1]:
              continue
            pragmaAnnotesDict[x][-1].append(prag)
          if prag not in regionDict.keys() or (prag in regionDict.keys() and regionDict[prag][-1] != 1):
              if res.qualifiers.get("sequence-region") and res.type in metaTypes:
                fields = res.qualifiers["sequence-region"][0].split(" ")
                regStr = writeMetaQuals({"sequence-region": res.qualifiers["sequence-region"][0]})
                regStr = regStr[regStr.find(" ") + 1:-1]
                fields = regStr.split(" ")
                regionDict[prag] = (int(fields[1]) - 1, int(fields[2]), annoteBit)
              elif res.qualifiers.get("is_circular") == ['True']:
                regionDict[prag] = (0, int(res.location.end), -1)
              else:
                regionDict[prag] = (0, int(res.location.end), 0)
          if suppressMeta <= 1 and res.type in metaTypes:
            qualsToAnnotes(metaAnnotesDict, res, prag)
            if suppressMeta == 1:
              orgDict[prag] = []
          if "Parent" in res.qualifiers.keys():
              seekParentDict[prag].append(0)  
        else: # Check if it's possible to resolve as a CompoundLocation feature
          if suppressMeta <= 1 and res.type in metaTypes:
            qualsToAnnotes(metaAnnotesDict, res, prag)
            if suppressMeta == 1:
              break  
          incInd = True
          if regionDict[prag][2] < 1:
            if res.qualifiers.get("is_circular") == ['True'] and int(res.location.end) > regionDict[prag][1]:
                regionDict[prag] = (0, int(res.location.end), -1)
            elif int(res.location.end) > regionDict[prag][1]: # Can't just max() or else we'll maybe overwrite -1 status
                regionDict[prag] = (0, int(res.location.end), 0)
          if res.id:
            for x in range(0, len(orgDict[prag])):
              if res.id == orgDict[prag][x].id:
                if orgDict[prag][x].type != res.type:
                  errOut += ("Line %d: Duplicate IDs in file but differing types. Cannot assume CompoundFeature/ join location, please resolve type descrepancy or de-duplicate ID %d.\n" % (lineInd, res.id))
                orgDict[prag][x].location = orgDict[prag][x].location + res.location
                incInd = False 
                break
          # If incInd is still true, then it's a unique feature, append to list
          if incInd:
            orgDict[prag].append(res)
            if "Parent" in res.qualifiers.keys():
              seekParentDict[prag].append(len(orgDict[prag])-1)
    
    orgDict, resolveErr = resolveParent(orgDict, seekParentDict)
    if resolveErr:
      errOut += resolveErr
    finalOrg = rAddDict(finalOrg, orgDict)

    # All features and pragmas should be read in by now, resolve any outstanding 
    # annotation or sequence associations

    for x in regionDict.keys():
      if x in seqDict.keys() and seqDict[x] != "": ## If x not in SeqDict, then a sequence-region pragma was made for organism with no features
        regionDict[x] = (0, len(seqDict[x]), 1)  # Make FASTA the final arbiter of region if present
    for x in regionDict.keys():
      if regionDict[x][2] == -1:
        errOut += "Organism %s: No sequence-region specified and last feature is labeled circular, unable to infer organism length.\n" % (x)       

    for x in finalOrg.keys():
      for i in finalOrg[x]:
        circ = False
        badIDs = []
        checkList = [i]
        for j in checkList: # make is_circular retroactively applicable to all features in a tree
          for k in j.sub_features:
            checkList.append(k)
          if j.qualifiers.get("is_circular") == ['True']:
            circ = True
            break
          if int(j.location.start) < regionDict[x][0] or int(j.location.end) > regionDict[x][1]:
            badIDs.append(j.id)
      if badIDs != [] and circ == False:
        errOut += "Organism %s: The following features fall outside of the specified sequence region: %s.\n" % (x, str(badIDs)[1:-1])     
         
    # By this point, all features and pragmas should be processed and resolved
    if errOut:
      outStream.write(errOut + "\n")
      raise Exception("Failed GFF Feature Parsing, error log output to stderr\n")

    # Construct a SeqRecord from all processed OrgIDs
    res = []
    for x in finalOrg.keys():
      finalOrgHeirarchy = []
      annoteDict = {}
      for pragma in pragmaAnnotesDict.keys():
        for vals in pragmaAnnotesDict[pragma]:
          if x in vals[1:]:
            annoteDict[pragma]=vals[0]
            break
      for pragma in metaAnnotesDict.keys():
        for vals in metaAnnotesDict[pragma]:
          if x in vals[1:]:
            if pragma in annoteDict.keys():
              if pragmaOverridePriority == False and pragBit < annoteBit:
                annoteDict[pragma]=vals[0]
                break
            else:
              annoteDict[pragma]=vals[0]
      for i in finalOrg[x]:
        if "Parent" not in i.qualifiers.keys():
          finalOrgHeirarchy.append(i)
          if i.type in metaTypes:
            if pragmaOverridePriority == False:
              for key in annoteDict.keys():
                if key not in finalOrgHeirarchy[-1].qualifiers.keys():
                  finalOrgHeirarchy[-1].qualifiers[key] = [annoteDict[key]]
            else:
              for key in annoteDict.keys():
                finalOrgHeirarchy[-1].qualifiers[key] = [annoteDict[key]]
            
      if seqDict[x]:
        if x in regionDict.keys():
          annoteDict["sequence-region"] = "%s %s %s" % (x, regionDict[x][0] + 1, regionDict[x][1])
          if len(seqDict[x]) < regionDict[x][1] - regionDict[x][0]:
            seqDict[x] += "?" * (regionDict[x][1] - regionDict[x][0] - len(seqDict[x]))
          else:
            seqDict[x] = seqDict[x][regionDict[x][0]:regionDict[x][1]]
        else:
          annoteDict["sequence-region"] = "%s 1 %s" % (x, int(len(seqDict[x])))
        seqDict[x] = str(seqDict[x])
      elif x in regionDict.keys():
        annoteDict["sequence-region"] = "%s %s %s" % (x, regionDict[x][0] + 1, regionDict[x][1])
        seqDict[x] = str(UnknownSeq(regionDict[x][1] - regionDict[x][0]))
      else: # Should actually no longer be reachable
        seqDict[x] = ""
      
      res.append(SeqRecord(Seq(seqDict[x]), x, "<unknown name>", "<unknown description>", None, finalOrgHeirarchy, annoteDict, None))
  
    return res

def printFeatLine(inFeat, orgName, source = 'feature', score = None, phase = None, outStream = sys.stdout, parents = None, codingTypes = ["CDS"]):
    for loc in inFeat.location.parts:
      line = orgName + "\t"
      if source:
        line += source + "\t"
      else:
        line += ".\t" 
      line += inFeat.type + "\t"
      startStr = str(min(loc.start, loc.end) + 1)
      endStr = str(max(loc.start, loc.end))
      if startStr[0] == "<":
        startStr = startStr[1:]
      if endStr[0] == ">":
        endStr = endStr[1:]
      line += startStr + "\t" + endStr + "\t"
      if score:
        line += str(score) + "\t"
      else:
        line += ".\t"
      if inFeat.location.strand == 0:
        line += ".\t"
      elif inFeat.location.strand == 1:
        line += "+\t"
      elif inFeat.location.strand == -1:
        line += "-\t"
      else: 
        line += "?\t"
      if inFeat.type in codingTypes: 
        if phase or phase == 0:
          line += str(phase) + "\t"
        else:
          line += "0\t"
      elif phase != 0:
        line += str(phase) + "\t"
      else:
        line += ".\t"
      if parents and "Parent" not in inFeat.qualifiers.keys():
        inFeat.qualifiers["Parent"] = parents.qualifiers["ID"]
      for qual in inFeat.qualifiers.keys():
        for keyChar in str(qual):
          if keyChar in "%,=;":
            encoded = str(hex(ord(keyChar)))
            line += "%" + encoded[2:].upper()
          else:
            line += keyChar
        line += "="
        if type(inFeat.qualifiers[qual]) != list:
          inFeat.qualifiers[qual] = [inFeat.qualifiers[qual]]
        for ind in range(0, len(inFeat.qualifiers[qual])):
          for valChar in str(inFeat.qualifiers[qual][ind]):
            if valChar in "%,=;":
              encoded = str(hex(ord(valChar)))
              line += "%" + encoded[2:].upper()
            else:
              line += valChar        
          if ind < len(inFeat.qualifiers[qual]) - 1:
            line += ","
          else:
            line += ";"
      
      outStream.write(line + "\n")  
    if type(inFeat) == gffSeqFeature and inFeat.sub_features: 
      for x in inFeat.sub_features:
        printFeatLine(x, orgName, x.source, x.score, x.phase, outStream, inFeat)

def gffWrite(inRec, outStream = sys.stdout, suppressMeta = 1, suppressFasta=True, codingTypes = ["CDS"], metaTypes = ["remark"], validPragmas = None, recPriority = True, createMetaFeat=None):

    writeFasta = False
    verOut = "3"
    firstRec = True

    if not inRec:
      outStream.write("##gff-version 3\n")
      return
    if type(inRec) != list:
      inRec = [inRec]
    for rec in inRec:
      if not isinstance(rec.seq, UnknownSeq):
        writeFasta = True
      
      seenList = []
      outList = {}
      if suppressMeta < 2:
        outList = rec.annotations
        #outStr = writeMetaQuals(rec.annotations) 
        #outStream.write(outStr)
        metaFeats = []
        for feat in rec.features:
          if feat.type in metaTypes:
            metaFeats.append(feat)
        for feat in metaFeats:
          for x in feat.qualifiers.keys():
            if recPriority == False or x not in outList.keys():
              outList[x] = " ".join(feat.qualifiers[x])
        
        if "gff-version" in outList.keys() and outList["gff-version"] != verOut:
          verOut = outList["gff-version"]
          outStream.write("##gff-version %s\n" % verOut)
          firstRec = False
        elif firstRec:
          outStream.write("##gff-version %s\n" % verOut)
          firstRec = False
        if validPragmas == None:
          outStr = writeMetaQuals(outList)
        else:
          whiteList = {}
          for x in outList.keys():
            if x in validPragmas:
              whiteList[x] = outList[x]
          outStr = writeMetaQuals(whiteList)
          outList = whiteList
        
        outStream.write(outStr)

      elif firstRec:
        outStream.write("##gff-version 3\n")

      foundMeta = False
      if createMetaFeat != None:
        for key in outList.keys():  # Change to GFF format qualifier dict
          outList[key] = [outList[key]]
        for feat in rec.features:
          if feat.type in metaTypes:
            foundMeta = True
            for key in outList.keys():
              if recPriority or key not in feat.qualifiers.keys():
                feat.qualifiers[key] = outList[key]
            break
        
        if not foundMeta:
          tempSeq = gffSeqFeature(FeatureLocation(0, len(rec.seq), 0), createMetaFeat, '', 0, 0, outList, None, None, None, '.', '.', "CPT_GFFParse") 
          printFeatLine(tempSeq, rec.id, source = tempSeq.source, score = tempSeq.score, phase = tempSeq.phase, outStream = outStream)

      for feat in rec.features:
          if suppressMeta > 0 and feat.type in metaTypes:
            continue  
          printFeatLine(feat, rec.id, source = feat.source, score = feat.score, phase = feat.phase, outStream = outStream)   
      firstRec = False 
    if writeFasta and not suppressFasta:
      outStream.write("##FASTA\n")
      for rec in inRec:
        rec.description = ""
        if not isinstance(rec.seq, UnknownSeq):
          SeqIO.write(rec, outStream, "fasta")
