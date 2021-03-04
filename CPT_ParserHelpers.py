# Copyright 2020-2021, Anthony Criscione
# Developed for the Center for Phage Technology, Texas A&M University
#
# Distributed under the BSD 3-Clause License, see included LICENSE file

# A collection of helper functions for use with the main GFF IO functions

from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation, CompoundLocation
from Bio.Seq import Seq, UnknownSeq
from collections import OrderedDict
try:
  from collections.abc import Iterable
except:
  from collections import Iterable
from gffSeqFeature import *

import sys

disallowArray = ["&", ",", ";", "="]
validArray = ["%26", "%2C", "%3B", "%3D"]
encoders = "ABCDEF1234567890"

validID = '.:^*$@!+_?-|'

def writeMetaQuals(qualList): 
    outLines = ""
    for x in qualList.keys():
      if x == "sequence-region":
        try:
          if isinstance(qualList[x], str):
            if qualList[x][0] == "(" and qualList[x][-1] == ")":
              fields = (qualList[x][1:-1]).split(" ") 
            else:
              fields = qualList[x].split(" ")
            if len(fields[0]) > 2 and fields[0][0] in ["'", '"'] and fields[0][0] == fields[0][-1]:
              fields[0] = fields[0][1:-1]
          
            if "%" in fields[1]:
              fields[1] = int(fields[1][:fields[1].find("%")])
            elif "," in fields[1]:
              fields[1] = int(fields[1][:fields[1].find("%")])
            else:
              fields[1] = int(fields[1])

            if "%" in fields[2]:
              fields[2] = int(fields[2][:fields[2].find("%")])
            else:
              fields[2] = int(fields[2])

          else:
            fields = qualList[x]
          outLines += "##sequence-region %s %d %d\n" % (fields[0], fields[1], fields[2])
        except:
          sys.stderr.write("Annotation Error: Unable to parse sequence-region in metadata feature. Value was %s" % (qualList[x]))
        
      elif x != "gff-version":
        outLines += "##%s" % (x)
        if isinstance(qualList[x], str):
          outLines += " %s" % (qualList[x].replace("\n", " "))
        elif isinstance(qualList[x], Iterable):
          for i in qualList[x]:
            outLines += " %s" % (str(i).replace("\n", " "))
        else:
          outLines += " %s" % (str(qualList[x]).replace("\n", " "))
        outLines += "\n"
    return outLines  

def validateID(idIn):
    badChar = []
    for x in idIn:
      if (ord(x) > 47 and ord(x) < 58) or (ord(x) > 64 and ord(x) < 91) or (ord(x) > 96 and ord(x) < 123) or (x in validID):
        continue
      else:
        if not(x in badChar):
          badChar.append(x)
    return badChar 

def replaceBadChars(qualIn):
    newQual = ""
    for x in qualIn:
      goodVal = True
      for y in range(0, len(disallowArray)):
        if x == disallowArray[y]:
          goodVal = False
          newQual += validArray[y]
      if goodVal:
        newQual += x
    return newQual

def validateQual(qualIn):
    badChar = []
    for x in qualIn:
      if x in disallowArray:
        if not(x in badChar):
          badChar.append(x)
    return badChar 

def rAddDict(lDict, rDict):
    for x in rDict.keys():
      val = lDict.get(x, [])
      val += rDict[x]
      lDict[x] = val
    return lDict

def checkCycle(orgDict):
  badOrgs = {}
  for org in orgDict.keys():
    for feat in orgDict[org]:
      if foundID(feat, feat.id):
        if org in badOrgs.keys():
           badOrgs[org].append(feat.id)
        else:
           badOrgs[org] = [feat.id]
          
  return badOrgs

def resolveParent(orgDict, indexDict):
  errOut = ""
  for org in indexDict.keys():
    for ind in indexDict[org]:
      for x in orgDict[org][ind].qualifiers['Parent']:
        for y in orgDict[org]:
          found = False
          if "ID" in y.qualifiers.keys() and x in y.qualifiers["ID"]:
            y.sub_features.append(orgDict[org][ind])
            found = True
            break
        if not found:
          errOut += ("Organism %s: Unable to find parent %s of feature %s\n" % (org, x, orgDict[org][ind].id))
  cycles = checkCycle(orgDict)
  if cycles.keys() != []:
    for x in cycles.keys():
      errOut += ("Organism %s: Cycle/ loop of features found involving feature IDs %s.\n" % (x, str(cycles[x])[1:-1]))
  if errOut != "":
    return None, errOut
  return orgDict, None

def foundID(featIn, topID):
  if not len(featIn.sub_features):
    return False
  for x in featIn.sub_features:
    if x.id == topID:
      return True
    for y in x.sub_features:
      if foundID(y, topID):
        return True
  return False

# A check for if an unencoded semicolon made it into the body of a qualifier value
# Sometimes occurs from manually edited Notes qualifiers
def encodeFromLookahead(remLine):
    for x in remLine:
      if x == "=":
        return False
      if x in ";,":
        return True
    return True # x == newline or EOF

def isNum(evalString):
  for x in range(0, len(evalString)): 
    if not(ord(evalString[x]) > 47 and ord(evalString[x]) < 58):
      return False
  return True

def qualsToAnnotes(inDict, feat, orgID):
  for x in feat.qualifiers.keys():
      if x not in inDict.keys():
        dictVal = " ".join(feat.qualifiers[x])
        outStr = writeMetaQuals({x: dictVal})
        if outStr == "":
          if x == "gff-version":
            outStr = feat.qualifiers[x][0]
          else:
            continue  
        else:
          outStr = outStr[outStr.find(" ") + 1:-1]
        inDict[x] = [[outStr, orgID]]
      else:
        contains = False
        for pragma in inDict.keys():
          for val in inDict[pragma]:
            if orgID in val[1:]:
              contains = True
              break 
        if not contains:
          dictVal = " ".join(feat.qualifiers[x])
          outStr = writeMetaQuals({x: dictVal})
          if outStr == "":
            if x == "gff-version":
              outStr = feat.qualifiers[x][0]
            else:
              continue 
          else:
            outStr = outStr[outStr.find(" ") + 1:-1]
          inDict[x].append([outStr, orgID])
  return inDict  

