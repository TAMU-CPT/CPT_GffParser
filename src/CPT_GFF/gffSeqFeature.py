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
        if qualifiers is None:
            try:
              qualifiers = OrderedDict()
            except:
              qualifiers = {}
        self.qualifiers = qualifiers
        self._id = id
        if "ID" in self.qualifiers.keys():
          self._id = self.qualifiers["ID"][0]
        elif id != "<unknown id>":
          self.qualifiers["ID"] = [id]
        if sub_features is None:
            sub_features = []
        self._sub_features = sub_features
        if ref is not None:
            # TODO - Deprecation warning
            self.ref = ref
        if ref_db is not None:
            # TODO - Deprecation warning
            self.ref_db = ref_db

    def _set_id(self, value):
        # TODO - Add a deprecation warning that the seq should be write only?
        self._id = value
        self.qualifiers["ID"] = [value]

    id = property(
        fget=lambda self: self._id,
        fset=_set_id,
        doc="The ID property, syncs with the qualifier field.",
    )

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
    noIDList = []
    lastCount = 0
    expectedParents = 0
    maxLoc = 0
    for feat in rec.features:
      if "Parent" in feat.qualifiers.keys():
        childList.append((convertSeqFeat(feat, defaultSource), [])) # Possible to have more than one parent
        expectedParents += len(feat.qualifiers["Parent"])
        #lastCount += childList[-1][1]
      elif feat.id and feat.id != "<unknown id>": # Do not accept the default value
        topList.append(convertSeqFeat(feat, defaultSource))
      else: 
        noIDList.append()
      maxLoc = max(maxLoc, feat.location.end)
    if deriveSeqRegion:
      rec.annotations["sequence-region"] = "%s 1 %s" % (rec.id, str(maxLoc))
    

    noEdit = False
    foundParCount = 0
    while not noEdit: 
      noEdit = True
      for childInd in range(0, len(childList)):
        for i in childList[childInd][0].qualifiers["Parent"]:
          nextChild = False
          for topFeat in topList:
            checkTree = [topFeat]
            for parCand in checkTree:
              nextPar = False
              checkTree += parCand.sub_features
              for foundPrior in childList[childInd][1]:
                if parCand.id == foundPrior:
                  nextPar = True
                  break
              if nextPar:
                break
              if i == parCand.id:
                parCand.sub_features.append(childList[childInd][0])
                childList[childInd]=(childList[childInd][0], childList[childInd][1] + [i])
                noEdit = False
                nextChild = True
                foundParCount += 1
                break
            if nextChild:
              break
      if noEdit and foundParCount < expectedParents:  
        badFeats = ""
        for x in childList:
          if len(x[0].qualifiers["Parent"]) != len(x[1]):
            badFeats += x.id + ", "
        sys.stderr.write("Unable to convert SeqRecord %s: could not find parents for features [%s]\n" % (rec.id, badFeats))

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
      
     
      
        
                 
     
