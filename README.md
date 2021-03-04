# CPT_GFFParser
An extensively featured Python parser for GFF format reading and writing

### Table of Contents  
- [Need](#Need)  
- [Goals](#Goals)  
- [Release](#Release)  
- [Library Contents](#library-contents)  
  - [CPT_GFFParser.gffParse](#CPT_GFFParsergffParse)
  - [CPT_GFFParser.gffWrite](#CPT_GFFParsergffWrite)
  - [gffSeqFeature.gffSeqFeature](#gffseqfeaturegffseqfeature)
  - [gffSeqFeature.convertSeqRec](#gffseqfeatureconvertSeqRec)
  - [gffSeqFeature.convertSeqFeat](#gffseqfeatureconvertSeqFeat)
- [Credit](#Credit)  
- [Contact Us](#contact-us)  
    

## Need
The Center for Phage Technology uses the Generic Feature Format extensively in our analysis, and Biopython is the most extensively integrated library within our own codebase. Unfortunately, the GFF format is not directly supported by Biopython, and as development of that project continued, we were no longer able to find an acceptable community solution to our need to support the GFF format. We encountered 2 main obstacles to implementing outside solutions: 
1. Biopython compatibility
   - It was imperative that any solution be able to interact cleanly with Biopython, as requiring extensive work to be done on the input would require rewrites across most of our codebase. However, Biopython has made significant changes over the years to its structure, and community solutions unfortunately had not kept pace. Additionally, as legacy versions of Biopython were cycled out of distribution, we could no longer rely on simply not upgrading as an ad-hoc solution,as losing our local copies of the packages could result in disaster.
2. Galaxy Compatibility
   - Most of the CPT's work is done on our instance of the Galaxy Project environment. Our Galaxy install distributes resources across local area machines and within one-off, limited permission environments. While there were more up-to-date GFF solutions out there, they also relied on a more concrete assumption of permanence and less finnicky writing permissions. We required a highly portable solution that was ambivalent about its working environment.

## Goals
Given a lack of available solutions, the decision was made to develop our own GFF library for parsing and writing out the format. 
1. Be a portable, light-weight, memory-only solution to reading in and writing out GFF format files. 
   - A Python library was a proven solution in the context of our other work.
2. Interface directly with SeqRecord and other Biopython structures to require minimal work on our main scripting files to be compatible.
   - Our prior solution had been BCBio's GFF functions. These generally used arrays of SeqRecord objects as their input arguments and output.
   - With this in mind, upon release, the overwhelming number of our tools were able to be made compatible simply by changing the import line, and replacing previous function names with gffParse and gffWrite. After tests on an intial batch of tools, the rest were successfully made compatible with nothing more than bash commands.
3. Related to the above, clean replacement of BCBio was desired, with minimal need to change assumptions about input arguments and especially without a need to change the XML structure of our existing Galaxy tools.
   - For the parser, this was done by maintaining that the first input was the GFF file handle, and the second, optional input was a dictionary object indexed on Organism ID that had SeqRecord objects as their value.
   - Similarly, the writer accepts the Biopython objects first, and then a writable handle or stream.
   - However, it should be noted that this is not an extension of BCBio; this code is not derived from that project, and a deliberate effort was made to not inspect their repository so as to avoid any influence. This is made without any claim either way as to quality, but our needs require an actively developed solution. We are grateful to Blue Collar Bioinformatics for their work, which sustained us all the way to 2020. 
4. More cleanly integrate aspects of the GFF format into existing Biopython based solutions.
   - The SeqFeature object and related constructs are very useful programming tools. Unfortunately, one major recent change to the library was to deprecate the sub_features functionality, in favor of using Compound Locations. 
   - While this reasoning is sound for formats Biopython directly supports, such as Genbank, it was less applicable to GFF3, which has heirarchy extremely baked into the concept, and more importantly for our needs, prior code has assumed a heirarchichal relationship between related features, of a wide array of types. 
   - Therefore, we decided to create a subclass of SeqFeature, gffSeqFeature, which reimplemented sub_features and could still be passed along as an array to unmodified SeqRecord objects. Additionally, guaranteed-present attributes, such as phase and score (Which are dedicated columns in the GFF spec) were added as attributes, rather than relying on equivalent qualifiers, which may be lost in the course of other bioinformatic work.
5. Finally, create a parser that is robust and informative at diagnosing failures and informing the user as to what went wrong.
   - A frequent demand from our users is validation scripts, to be run on datasets to ensure they conform to some kind of spec before doing extensive work and possibly wasting that time if the initial data is poorly structured.
   - We desired to extend this experience to the parser, and create a tool that, as best as possible, could withstand reading in improperly formatted GFF files until the EOF, and output a comprehensive list of problematic lines and their errors. It hoped that the errorlog this program creates will have enough information that the user will only need to go back for one round of changes.

## Release
gffParse and gffWrite were initially released at the end of Summer, 2020 on the CPT Galaxy server, effective across all of our tools interfacing with GFF files. Thanks to the testing and patience of our users and students, it was a mostly stable library by the end of September, and almost fully featured by that November. February 2021 marks the most recent round of feature implementation, primarily focused on metadata features and metadata directives the GFF spec names pragmas (The initial release supported the ##FASTA pragma, but the others were not deemed a priority at the time). Additionally, in anticipation of a public standalone release, more arguments were added to the primary IO functions in order to help users fit the tool to their own needs. For example, the CPT frequently finds it highly useful in our workflow to strip as much metadata as possible, however we have now implemented a variety of options for writing out and dynamically creating metadata features and annotations if desired.

## Library Contents
### CPT_GFFParser.gffParse
The primary function for reading in GFF files. Will return a list of SeqRecord objects, with gffSeqFeature objects as their .feature lists.

`gffParse(gff3In, base_dict = {}, outStream = sys.stderr, codingTypes=["CDS"], metaTypes = ["remark"], suppressMeta = 2, pragmaPriority = True, pragmaOverridePriority = True):`
- gff3In --- source file handle, will accept any object that implements a .readlines() method or other list-like access to file lines
- base_dict --- Additional SeqRecord information. Keys are OrganismIDs and values are SeqRecords. For BCBio backwards compatibility.
- outStream --- output filestream or stringstream for the errorlog to be passed, if any parsing errors are encountered
- codingTypes --- list of feature types where a non-'.' phase value is expected, passed along to lineAnalysis
- metaTypes --- list of metadata feature types. Features of this type will be affected by the remaining arguments
- suppressMeta --- Suppress metadata fields. Integer value, where: 
  - 0 == no suppression, all metadata from features and pragmas will be read and output to the SeqRecord as .annotation entries.
  - 1 == As above, but metadata features will not be entered into the SeqRecord's feature list after their metadata is recorded and entered into the SeqRecord.annotation
  - 2 == Total suppression, no metadata features will even be processed, and no pragmas except those related to sequence length (##FASTA and ##sequence-region) and ##gff-version (required by GFF spec) will be utilized.
- pragmaPriority --- In cases where pragmas and metadata features disagree/conflict, pragmas will take precedence for creating SeqRecord.annotation value if true, else the feature will.
- pragmaOverridePriority --- Similar to above, in the event of a conflict between metadata features and pragmas, the pragma's value will override the metadata gffSeqFeature.qualifier value with the pragma's own. This will force the metadata and pragmas to sync, and avoid future discrepancies. Should only be used with pragmaPrority also set to True

### CPT_GFFParser.gffWrite
The primary function for writing out GFF files. 

`gffWrite(inRec, outStream = sys.stdout, suppressMeta = 1, suppressFasta=True, codingTypes = ["CDS"], metaTypes = ["remark"], validPragmas = None, recPriority = True, createMetaFeat=None)`
- inRec --- The input, can either be one SeqRecord object or a list of them. Expects the features to be of gffSeqFeature type, if working with a vanilla Biopython record (For example, after having used SeqIO to read in a Genbank file), please see `gffSeqFeature.convertSeqRec` below.
- outStream --- The output location, can be a file handle, stringstream, or anything else with a .write method implemented.
- suppressMeta --- Suppress metadata fields. Integer value, where: 
  - 0 == no suppression, all metadata from features of a type in metaTypes and annotations from SeqRecord.annotation will be read and output to the GFF as ##pragma entries, where the key will be the pragma and the values will be joined by " " and output on the same line (newlines replaced by " ").
  - 1 == As above, but metadata features will not be written out to file after their annotation information is retrieved.
  - 2 == Total suppression, no metadata features or SeqRecord.annotations will even be processed, and no pragmas or metadata features will be in the final output. The two exceptions are at least one ##gff-version pragma will be created at the first line, as per the GFF spec requirements, and if a FASTA is supplied and suppressFasta is set to False, a ##FASTA directive and corresponding sequence(s) will be created at the end of the file.
- suppressFasta --- If True, do not write the Fasta sequences of this SeqRecord list at the end of the GFF
- codingTypes --- A list of feature types where an integer is expected in the phase field, even if that integer would be the default value of 0. Features of types not in this list will have their phase written as the default "." empty value IF their .phase is 0.
- metaTypes --- A list of feature types that are metadata information. Features of these types will have their qualifiers counted as annotations, and they will be affected by the actions of suppressMeta above. 
- validPragmas --- A whitelist of pragmas to output, in cases where some, but not all, annotations would be desired. Setting to None will allow allow all annotations. Setting to empty list will allow no annotations (Except as discussed above in suppressMeta, where #gff-version and ##FASTA are compulsory). Note that this affects pragmas only, and metadata *features* will not be affected by a choice here.
- recPriority --- In cases where a SeqRecord.annotation and a qualifier in a metadata feature conflict, the SeqRecord will take priority if set to true, otherwise the metadata feature will.
- creatMetaFeat --- A string input, where if it is set to None then no metadata feature will be created, otherwise a metadata feature of that type will be created, and all annotations wil be written out as its qualifiers. For example, createMetaFeat="remark" will cause a feature of type remark to be created with appropriate FeatureLocation and qualifiers. If there already exists a feature of that type in the output, then that feature will simply be updated.

### gffSeqFeature.gffSeqFeature
A subclass of SeqFeature. Minimal editing has been done to it, so it should still be largely compatible with any operations that are done on a regular SeqFeature. However the sub_features property has been reimplemented, so care must be taken that analysis done on "top-level" features cascades down to sub-features, if need be. See the documentation on Biopython's SeqFeature for most needs, only the differences are discussed here:

`gffSeqFeature(SeqFeature.SeqFeature):
    def __init__(self, location=None, type="", location_operator="", strand=None, id="<unknown id>", qualifiers=None, sub_features=None, ref=None, ref_db=None, phase=0, score=0.0, source="feature"):`
Three new fields have been appended, phase, score, and source, columns 8, 6, and 2 respectively in GFF output. All three can be accessed as .properties of a gffSeqFeature object. Additionally, sub_features will now accept a list of gffSeqFeatures as input, rather than informing the user this feature has been deprecated.

`    def _shift(self, offset):`
Behaves identically to \_shift in a SeqFeature object, but reimplemented because the output has to be explicitly cast as a gffSeqFeature.

`    def translate(self, parent_sequence, table="Standard", start_offset=None, stop_symbol="*", to_stop=False, cds=None, gap=None):`
Behaves identically to regular SeqFeature, except the .phase property is checked along with seeing if there's a codon_start qualifier. start_offset will still override both, and codon_start will take precedence over .phase

### gffSeqFeature.convertSeqRec
A function for converting the members of a SeqRecord's SeqFeature list into gffSeqFeature objects. Will output a list of SeqRecord objects. If your input file was not a GFF then you probably want to run this before sending your SeqRecord to gffWrite. Will attempt to pair "source", "codon_start", "score" and "Parent" qualifiers in a SeqFeature with the appropriate gffSeqFeature property.

As long your SeqFeature objects have "Parent" qualifiers whose value corresponds to the .id property of another SeqFeature in the list, then it will correctly construct a heirarchy. Note that this script will fail if you construct a series of features that parent each other in a loop or cycle, or if you create a feature who lists itself as a Parent (self-loop). These qualifiers are allowed in Biopython because there is no formal heirarchy, but GFF must enforce it.

`convertSeqRec(inRec, defaultSource = "gffSeqFeature", deriveSeqRegion = True, createMetaFeat = None):`
- inRec --- A seqRecord or list of SeqRecords to have their features converted.
- defaultSource --- A string input for the source field of the gffSeqFeature (column 2 in GFF output). If a "source" qualifier is present in the feature, that will be used instead.
- deriveSeqRegion --- If True, the script will derive a ##sequence-region pragma/annotation based off of the "largest" FeatureLocation.end encountered while reading the features in. It will then add this key/value to the SeqRecord.annotation dictionary. 
- createMetaFeat --- A string input, where if it is set to None then no metadata feature will be created, otherwise a metadata feature of that type will be created, and all annotations wil be written out as its qualifiers. For example, createMetaFeat="remark" will cause a feature of type remark to be created with appropriate FeatureLocation and qualifiers. If there already exists a feature of that type in the output, then that feature will simply be updated.

### gffSeqFeature.convertSeqFeat
Primarily used by gffSeqFeature.convertSeqRec, but can be accessed as a standalone function if there is a need to convert only individual features detached from a SeqRecord object.

`convertSeqFeat(inFeat, defaultSource = "gffSeqFeature"):`
- inFeat --- A single FeatureLocation object
- defaultSource --- If there is no "source" qualifier in this feature's qualifiers, this string will be set as the source (column 2 in GFF output)

## Credits
Center for Phage Technology,
Anthony Criscione, 2020-2021

Distributed under the BSD 3-Clause license, included in the LICENSE file of this same directory

gffSeqFeature developed in accordance with the Biopython License Agreement, and similarly licensed under the BSD 3-Clause License

## Contact Us
For issues directly related to coding errors or other library problems, we encourage you to submit an issue via the issues tab above. Please provide relevant information such as error messages, function arguments, and input data.

For other inquiries, contact info for the Center for Phage Technology can be found at cpt.tamu.edu
