# CPT_GffParser
An extensively featured Python parser for GFF format reading and writing

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
GFFParse and GFFWrite were initially released at the end of Summer, 2020 on the CPT Galaxy server, effective across all of our tools interfacing with GFF files. Thanks to the testing and patience of our users and students, it was a mostly stable library by the end of September, and almost fully featured by that November. February 2021 marks the most recent round of feature implementation, primarily focused on metadata features and metadata directives the GFF spec names pragmas (The initial release supported the ##FASTA pragma, but the others were not deemed a priority at the time). Additionally, in anticipation of a public standalone release, more arguments were added to the primary IO functions in order to help users fit the tool to their own needs. For example, the CPT frequently finds it highly useful in our workflow to strip as much metadata as possible, however we have now implemented a variety of options for writing out and dynamically creating metadata features and annotations if desired.

## Library Contents
### CPT_GFFParser.GFFParse

### CPT_GFFParser.GFFWrite

### gffSeqFeature.gffSeqFeature

### gffSeqFeature.convertSeqRec

### gffSeqFeature.convertSeqFeat

## Credits
Center for Phage Technology,
Anthony Criscione, 2020-2021
