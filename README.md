Differential Expression
=======================
These Python, R and BASH scripts have been used to process BAM files containing Aedes Aegypti reads aligned to a reference genome obtained from VectorBase. Many of the scripts in this repository cannot be properly used on another computer, since many of BAM files are too large to upload.

These BAM files each contain a set of "reads" which are digitized DNA sequences output from a DNA sequencer. The sequences correspond to the set of mRNA molecules extracted from an organism under a particular experimental condition. The sequences in the BAM files have been mapped to an "reference genome" and thus contain the corresponding numeric positions (also called "chromosomal coordinates") to which they map in the reference.

copy_bam.sh: moves the BAM files stored on the local computer to an external hard drive labeled JAXI. The files are then processed with a program called featureCounts. featureCounts takes every read in each BAM file supplied to it and check whether that read overlaps with any annotated sequences in an "annotated genome" which contains several metadata for each known sequence in the organisms genome including the starting and ending chromosomal coordinates as well as its length.
