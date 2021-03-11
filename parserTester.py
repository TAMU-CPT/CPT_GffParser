#!/usr/bin/env python
import argparse
import logging
import sys
from CPT_GFFParser import gffParse, gffWrite
from gffSeqFeature import convertSeqRec
from Bio.Seq import Seq, UnknownSeq
from Bio import SeqIO


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract specified qualifers from features in GFF3', epilog="")
    parser.add_argument('gff3', type=argparse.FileType("r"), help='GFF3 File')
    args = parser.parse_args()
    recs = gffParse(args.gff3, suppressMeta = 0, pragmaPriority = True, pragmaOverridePriority = True)
    featList = []
    for x in recs[0].features:
      featList.append(x)
      for y in x.sub_features:
        featList.append(y)
        for z in y.sub_features:
          featList.append(z)

    for x in featList:
      x.sub_features = []
    recs[0].features = featList
    recs = convertSeqRec(recs, defaultSource = "ConvertSeqRec")

    gffWrite(recs, suppressMeta = 0, createMetaFeat="remark")
