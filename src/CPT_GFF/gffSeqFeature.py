# Copyright 2020-2021, Anthony Criscione
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
from Bio import SeqFeature
from Bio.SeqFeature import FeatureLocation, CompoundLocation
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq, UnknownSeq
from collections import OrderedDict
import sys

class gffSeqFeature(SeqFeature.SeqFeature):
    def __init__(
        self,
        location=None,
        type="",
        location_operator="",
        strand=None,
        id="<unknown id>",
        qualifiers=None,
        sub_features=None,
        ref=None,
        ref_db=None,
        phase=0,
        score=0.0,
        source="feature"
    ):
        """Reimplementation of SeqFeature for use with GFF3 Parsing
        Does not remove the sub_feature functionality, as unlike
        Genbank, this is baked into the core concept of GFF
        """
        if (
            location is not None
            and not isinstance(location, FeatureLocation)
            and not isinstance(location, CompoundLocation)
        ):
            raise TypeError(
                "FeatureLocation, CompoundLocation (or None) required for the location"
            )
        self.location = location
        self.type = type
        self.phase = phase
        self.score = score
        self.source = source
        if location_operator:
            # TODO - Deprecation warning
            self.location_operator = location_operator
        if strand is not None:
            # TODO - Deprecation warning
            self.strand = strand
        self.id = id
        if qualifiers is None:
            try:
              qualifiers = OrderedDict()
            except:
              qualifiers = {}
        self.qualifiers = qualifiers
        if sub_features is None:
            sub_features = []
        self._sub_features = sub_features
        if ref is not None:
            # TODO - Deprecation warning
            self.ref = ref
        if ref_db is not None:
            # TODO - Deprecation warning
            self.ref_db = ref_db

    def _get_subfeatures(self):
        """Get function for the sub_features property (PRIVATE)."""
        try:
            return self._sub_features
        except AttributeError:
            return None

    def _set_subfeatures(self, value):
        """Set function for the sub_features property (PRIVATE)."""
        if isinstance(value, list):
            self._sub_features = value
        else:
            raise ValueError("sub_feature must be a list of gffSeqFeature objects")

    sub_features = property(
        fget=_get_subfeatures,
        fset=_set_subfeatures,
        doc="Sub-features for GFF Heirarchy",
    )

    def _shift(self, offset):
        """Return a copy of the feature with its location shifted (PRIVATE).
        The annotation qaulifiers are copied.
        """
        for x in self.sub_features:
          x._shift(offset)  
        return gffSeqFeature(
            location=self.location._shift(offset),
            type=self.type,
            location_operator=self.location_operator,
            id=self.id,
            qualifiers=OrderedDict(self.qualifiers.items()),
            sub_features=self.sub_features,
            phase=self.phase,
            score=self.score,
            source=self.source
        )

    def translate(
        self,
        parent_sequence,
        table="Standard",
        start_offset=None,
        stop_symbol="*",
        to_stop=False,
        cds=None,
        gap=None,
    ):
        """
          Identical to the implementation found in 
          Biopython SeqFeature, but will use .phase value instead
          if start_offset is not set and start_codon is not present

          Deferred to codon_start under reasoning that some bioinformatic scripts
          may edit the codon_start field, but not change the .phase value
        """
        # see if this feature should be translated in a different
        # frame using the "codon_start" qualifier
        if start_offset is None:
            try:
                start_offset = int(self.qualifiers["codon_start"][0]) - 1
            except KeyError:
                start_offset = self.phase

        if start_offset not in [0, 1, 2]:
            raise ValueError("The start_offset must be 0, 1, or 2. The supplied value is '%s'. Check the value of either the codon_start qualifier, the .phase property, or the start_offset argument" % (start_offset))

        feat_seq = self.extract(parent_sequence)[start_offset:]
        codon_table = self.qualifiers.get("transl_table", [table])[0]

        if cds is None:
            cds = self.type == "CDS"

        return feat_seq.translate(
            table=codon_table,
            stop_symbol=stop_symbol,
            to_stop=to_stop,
            cds=cds,
            gap=gap,
        )

def convertSeqFeat(inFeat, defaultSource = "gffSeqFeature"):
  featLoc = inFeat.location
  IDName = inFeat.id
  qualDict = inFeat.qualifiers
  parentCands = inFeat.qualifiers.get("Parent", [])
  for x in parentCands:
    if x == inFeat.id: # Cannot allow self-loops
      raise Exception("Cannot convert SeqRecord, feature %s lists itself as a parent feature" % (cand.id))
  if "codon_start" in inFeat.qualifiers.keys():
    phaseIn = int(inFeat.qualifiers["codon_start"][0])
  else:
    phaseIn = 0
  if "score" in inFeat.qualifiers.keys():
    scoreIn = float(inFeat.qualifiers["score"][0])
  else:
    scoreIn = "."
  if "source" in inFeat.qualifiers.keys():
    sourceIn = inFeat.qualifiers["source"][0]
  else:
    sourceIn = defaultSource 

  return gffSeqFeature(featLoc, inFeat.type, '', featLoc.strand, IDName, qualDict, [], None, None, phaseIn, scoreIn, sourceIn)

def convertSeqRec(inRec, defaultSource = "gffSeqFeature", deriveSeqRegion = True, createMetaFeat = None):
  # Assumes an otherwise well-constructed SeqRecord that just wants to replace its features with gffSeqFeatures
  if not isinstance(inRec, list):
    inRec = [inRec]
  
  outRec = []
  for rec in inRec: 
    topList = []
    childList = []
    noPairList = [] # Features with no ID that shouldn't be analyzed for sub-feature relations
    lastCount = 0
    maxLoc = 0
    for feat in rec.features:
      if "Parent" in feat.qualifiers.keys():
        childList.append((convertSeqFeat(feat, defaultSource), len(feat.qualifiers["Parent"]))) # Possible to have more than one parent
        lastCount += childList[-1][1]
      elif feat.id and feat.id != "<unknown id>": # Do not accept the default value
        topList.append(convertSeqFeat(feat, defaultSource))
      else:
        noPairList.append(feat)
      maxLoc = max(maxLoc, feat.location.end)
    if deriveSeqRegion:
      rec.annotations["sequence-region"] = "%s 1 %s" % (rec.id, str(maxLoc))
    
    rerunList = True
    while rerunList:
      rebuildList = False
      for ind in range(0, len(childList)): # Check for subfeatures of subfeatures first
        child = childList[ind]
        foundItem = child[1]
        for cand in childList:
          if foundItem > 0:
            for childID in child[0].qualifiers["Parent"]:
              if cand[0].id == childID:
                cand[0].sub_features.append(child[0])
                foundItem -= 1
                childList[ind] = (child[0], foundItem)
                rerunList = True
          elif foundItem == 0:
            break
      
    lastCount = -1
    thisCount = -1
    while lastCount != 0: # This shouldn't need to actually loop
      thisCount = 0
      popList = []
      for child in childList: 
        foundItem = child[1]
        for cand in topList:
          if foundItem > 0:
            for childID in child[0].qualifiers["Parent"]:
              if cand.id == childID:
                cand.sub_features.append(child[0])
                foundItem -= 1
          elif foundItem == 0:
            break
        if foundItem > 0:
          popList.append((child[0], foundItem))
          thisCount += popList[-1][1]
      childList = popList
      if thisCount != 0:
        popList = []
        lastCount = thisCount
        thisCount = 0 
      elif thisCount == lastCount or thisCount > 0:
        badIDs = []
        for x in childList:
          badIDs.append(x[0].id)
        outStr = ", ".join(badIDs)
        sys.stderr.write("Unable to convert SeqRecord %s: could not find parents for features [%s]\n" % (rec.id, outStr))
        sys.stderr.write("Note that this error will also occur if sub_feature relationships between features ever form a cycle/loop.\n")
        raise Exception("Could not convert features of SeqRecord %s to gffSeqFeature format, see stderr\n" % (rec.id)) 
      else:
        break

    if createMetaFeat:
      qualDict = {}
      for x in res.annotations.keys():
        outVal = ""
        if isinstance(res.annotations[x], list):
          outVal = " ".join(res.annotations[x])
        else:
          outVal = str(res.annotations[x])
        outVal = outVal.replace("\n"," ")
        qualDict[x] = [outVal]
      topList.append(gffSeqFeature(FeatureLocation(0, maxLoc), createMetaFeat, '', 0, IDName, qualDict, [], None, None, 0, ".", defaultSource))
    topList = sorted(topList, key=lambda feature: feature.location.start)
    rec.features = topList
    outRec.append(rec)
  return outRec
      
     
      
        
                 
     
