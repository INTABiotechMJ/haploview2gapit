 The input file is an output from the haploview program, a comma-separated values (csv) file containing the cultivar name in the first column and in the consecutive column, an haplotype block (sequence with more than one base) or a marker (single base sequence). For each pair of culitvar / haplotype or marker there's a sequence or a base.

First, incomplete sequences are searched in every cell. Those which count of missing data (indicated with a "-") is more than 0 and less than 10% of the sequence length. Also the sequence must have at least one missing data position. Notice that single base markers do not fit this criteria.

It is intended to replace the missing data positions in the incomplete sequences with a nucleotide that is in the same position in the complete sequences. All incomplete sequences are aligned (using a global pairwise alignmen algorithm) against each one of the full sequences (those without missing data) of the same haplotype block. 

A cutoff value is calculated as the quantity of nucleotides in the sequence (with no missing data) multiplied by 0.8. 

For all the complete sequences alignments that match previous criteria, the nucleotides located in the same position as the first missing data in the incomplete sequence are extracted. If all correspond to the same DNA base, the missing data is replaced with that.

Once all incomplete sequences are scanned and no replacement is done, a file for input in Gapit is generated. The first column correspond to a cultivar name, and a new column is added for each haplotype and for each different sequence from that haplotype. A 1 is placed if the sequence of the cultivar matches the sequence of the haplotype in the column name, and a 0 otherwise.