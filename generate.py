# -*- coding: utf-8 -*-
#!/usr/bin/env python
'''
TODO: Description of the script
'''
import argparse
import pandas as pd
import csv
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from operator import itemgetter
from math import floor

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Input file in .csv",  required=True)
parser.add_argument("-o", "--output", help="Output file in .csv",  required=True)
parser.add_argument("-v","--verbosity", help="increase output verbosity")

args = parser.parse_args()
def make_log(output):
    if args.verbosity:
        print output

df = pd.read_csv(args.input,sep=",")

all_data = []
allacc=[]
has_missing_data = True
count = 0
while has_missing_data:
    has_missing_data = False
    #for each marker
    for column in df:
        if column == "Marker":
            continue
        make_log("*"*20)
        make_log("Marker: " + column)
        for k,sequence in df[column].iteritems():
            # we have a missing value
            if not isinstance(sequence, basestring):
                print 'Sequence ' + str(sequence) + ' is not a string'
                continue
            if sequence.count('-') * 100 / len(sequence) > 10:
                make_log('More than 10% missing spots, skipping')
                continue
            if len(sequence.replace("-","")) == 0:
                make_log("Empty sequence, skipping")
                continue
            if not "-" in sequence:
                make_log("Full sequence, skipping")
                continue
            make_log("-"*20)
            cutoff = floor(len(sequence.replace("-","")) * 0.8)
            make_log("Sequence: "+sequence+ " - cut-off value:" + str(cutoff) + " - comparing to:")
            # search for position
            position = sequence.find("-")
            consensus = {}
            for v2 in df.loc[: ,column]:
                if v2.find("-") > -1:
                    continue
                alignments = pairwise2.align.globalxx(sequence, v2)
                if alignments:
                    max_score = max(alignments,key=itemgetter(2))[2]
                    #make_log(v2+ " - Alignment score: " + str(max_score))
                    if max_score > cutoff:
                        nucleotide = v2[position:position+1]
                        make_log("Accepted, alt nucleotide: " + nucleotide)
                        if nucleotide in consensus:
                            consensus[nucleotide] += 1
                        else:
                            consensus[nucleotide] = 1
            #get the consensus
            if len(consensus) == 1:
                has_missing_data = True
                #max_nucleotide = max(consensus.iteritems(), key=operator.itemgetter(1))[0]
                unique_nucleotide = consensus.keys()[0]
                make_log("Found unique nucleotide: " + unique_nucleotide)
                df.loc[k,column] = sequence[:position] + unique_nucleotide + sequence[position+1:]
            elif len(consensus) == 0:
                make_log("No alternatives found")
            else:
                make_log("More than one nucleotide, skipping: " + str(consensus))
            #replace in current value

#generate tuples (haplo / sequence) for each sequence
for column in df:
    for k,v in df.iterrows():
        if column == "Marker":
            allacc.append(v.Marker)
            continue
        seq = v[column]
        res = (column, seq)
        if not res in all_data:
            all_data.append(res)

#res = {}
#for re in all_data:
#    mark, seq = re
#    for acc in df["Marker"]:
#        cseq = df[(df.Marker == acc)].loc[:,mark].iloc[0]
#        if cseq == seq:
#            res[(acc,mark,seq)] = cseq == seq
#for col in all_data:
#    mark, trash = re
#    for r in res:
#        pass
#
outcsv = open(args.output, 'w')
writer = csv.writer(outcsv, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
li = []
li.append("Marker")
for ad in all_data:
    marker, sequence = ad
    li.append(marker + "_" + sequence)
writer.writerow(li)
for acc in allacc:
    li = []
    li.append(acc)
    for ad in all_data:
        marker, sequence = ad
        cur_seq = df[(df.Marker == acc)].loc[:,marker].iloc[0]
        if cur_seq == sequence:
            li.append("1")
        else:
            li.append("0")
    writer.writerow(li)
