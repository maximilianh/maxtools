#!/usr/bin/python 

# input: bed file with exons as exported from ucsc
# output: bed file with exons numbers reversed if strand="-"

from sys import *

f = open(argv[1], "r")
maxExon = {}
for l in f:
    fs = l.split()
    fs2 = fs[3].split("_")
    gene,exon = fs2[0],fs2[2]
    
    maxExon[gene] = int(exon)

f = open(argv[1], "r")

for l in f:
    fs = l.split()
    gene = fs2[0]
    exon = fs2[2]
    fs2 = fs[3].split("_")
    if fs[5]=="-":
        name = gene + "_exon_" + str(maxExon[gene] - int(exon))
    else:
        name = gene + "_exon_" + exon
    print "\t".join([fs[0], fs[1], fs[2],name,fs[4],fs[5]])


