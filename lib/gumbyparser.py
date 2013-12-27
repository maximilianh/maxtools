#!/usr/bin/python

from sys import *
from re import * 
from optparse import OptionParser
import Fasta
import alignment

class gumbyBlock:
    # easier to handle as a class
    # baseseq is always first in seqs list
    def __init__(self,number, score, pval, seqs):
        self.number=number
        self.score = score
        self.pval=pval
        self.seqs=seqs
    def __repr__(self):
        lines = []
        lines.append(  " * gumbyResult %d" % self.number)
        lines.append("score %d, pval %e " % (self.score, self.pval))
        for s in self.seqs:
            lines.append( str(s))
        return "\n".join(lines)


# -------- FUNCTIONS ------------------
def procData(baseSeq, exons, no, seqDict, pos, pval, length, score):
    if len(seqDict)==0:
        return []

    if baseSeq not in seqDict:
        stderr.write("error: your baseseq name is not in gumby result. remember that gumby keeps only first word of seq name\n")
        sys.exit(1)
    print seqDict[baseSeq]
    if overlapped(pos[baseSeq], exons, baseSeq):
        stderr.write("warning: dropping complete block with sequence %s:%s because baseSeq has overlapping exon annotation.\n" % (baseSeq, pos[baseSeq]))
        return []
    if seqDict[baseSeq].nucl.count("-")==len(seqDict[baseSeq].nucl):
        stderr.write("warning: dropping complete block with sequence %s:%s because baseSeq contains only '-'-characters\n" % (baseSeq, pos[baseSeq]))
        return []
    if seqDict[baseSeq].nucl.count("N")==len(seqDict[baseSeq].nucl):
        stderr.write("warning: dropping complete block with sequence %s:%s because baseSeq contains only N-characters\n" % (baseSeq, pos[baseSeq]))
        return []
    seqs = []
    seqs.append(seqDict[baseSeq])

    for n,s in seqDict.iteritems():
        if n==baseSeq:
            continue
        seqs.append(s)

    gb =  gumbyBlock(no, score, pval, seqs) 
    return [gb]

def overlapped(pos, exons, baseSeq):
    f1name, f1start, f1end = pos
    if f1name != baseSeq:
        return False
    for e in exons:
        f2start, f2end = e 
        # print "checking %d -- %d, %d" % (start, f2start, f2end)
        result = (( f2start <= f1start and f2end > f1start) or \
            (f2start < f1end and f2end >= f1end) or (f2start >= f1start and f2end <= f1end))
        if result == True:
            return True
    return False

# ----- MAIN -------------------
def parseGumby(gumbyFile, exonFile, baseSeq):
# parses gumbyFile, removes things that overlap exons and gumbies that consist only of gaps on baseSeq
# returns a list of gumbyBlocks

    infile = open(gumbyFile, "r")

    exons = []
    if exonFile!=None:
        fh = open(exonFile, "r")
        for l in fh:
            fs = l.split()
            if fs[0].lower()!=baseSeq:
                continue
            exons.append([ int(fs[3]), int(fs[4]) ] )
    # print exons

    re1 = compile("[a-z]+[ ]+[0-9]+[ ]+[0-9]+")
    seqs = {}
    pos = {}
    i = -1

    resultLst = alignment.Alignment()
    for l in infile:
        l = l.strip()
        l = l.replace("*","-")
        l = l.replace("<", "-")
        l = l.replace(">", "-")
        if l.startswith("start"):
            if i!=-1:
                resultLst.extend(procData(baseSeq, exons, i, seqs, pos, pval, length, score))
            f = l.split()
            pval = float(f[-1])
            length = int(f[6].strip(","))
            score = int(f[8].strip(","))
            i+=1
            seqs={}

        if re1.match(l):
            f = l.split()
            name = f[0]
            start = int(f[1])-1
            end = int(f[2])-1

            seq = f[3]
            if name not in seqs:
                faseq = Fasta.FastaSeq(name, seq)
                faseq.chrom = name
                faseq.start = start
                faseq.end = end
                seqs[name] = faseq
            else:
                faseq = seqs[f[0]] 
                faseq.nucl += f[3]
            pos[name] = (name, start,end)

    resultLst.extend(procData(baseSeq, exons, i, seqs, pos, pval, length, score))
    return resultLst

# ---- DEBUGGING -----
 
#blocks = parseGumby("test/gumbyparser.gumby", "test/gumbyparser.bed", "oryzias")
#for b in blocks:
    #print b
