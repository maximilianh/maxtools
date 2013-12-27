from array import *
import sys
import copy

""" put bed files onto lines = dict with y-deplacement => list of bed features """
def putOntoLines(beds, seqlen):
    lines = {}
    nextline = array('i', (seqlen+20) * [0])
    for h in beds :
            lastline = 0
            maxline = 0
            if h.start < seqlen:
                for i in range(h.start,max(h.end, h.start+len(h.name)+1)):
                    if i >= seqlen:
                        print "warning: cannot draw, feature too long:",h,"\n"
                        h.end = seqlen
                        break
                    line = nextline[i]
                    nextline[i] += 1 
                    maxline = max(maxline, line)
                lines.setdefault(maxline, [])
                lines[ maxline ].append(h)
    return lines

def ungappedToGapped(seqstr, beds, minlimit=-1, minend=-1):
    """ convert coords of beds and also minlimit and minend if (!= -1)"""
    """ return a three tuple of these """
    gapped = 0
    ungapped = 0
    newbeds = copy.copy(beds)
    maxpos = 1 + len(seqstr)
    uToG = array("L", maxpos * [0]) # ungapped to gapped array

    for c in seqstr:
        if c!="-":
            uToG[ungapped] = gapped
            ungapped += 1
        gapped += 1

    for b in newbeds:
        if b.start < maxpos and b.end < maxpos:
            b.start=uToG[b.start]
            b.end = uToG[b.end]
        else:
            sys.stderr.write( "warning: cannot put feature onto seq, feature pos exceeds alignment len:%s\n"% str(b))
    if minlimit==-1:
        minlimit=0
    if minend==-1:
        minend=0
    return (newbeds, uToG[minlimit], uToG[minend])

def paint(seqs, lines, minPos, maxPos):
    """ return ascii logo of seqs and features between min and max """
    seqlen = len(seqs[0].nucl)
    textLines = []
    # sequences
    header = ""
    for s in seqs:
        header = "%-25s : " % (s.id)
        textLines.append(header + s.nucl[minPos:maxPos])
    #print "head", len(header), "min", minPos, "max", maxPos
        

    # features
    for i in range(0, len(lines)):
        line = lines[i]
        size = maxPos-minPos+30
        linestr1=array("c", size * [" "])
        linestr2=array("c", size * [" "])
        for b in line:
            # get corrected positions 
            if b.end > maxPos:
                continue
            if b.start < minPos:
                continue
            start = b.start - minPos
            end = b.end - minPos

            # draw line
            for x in range(start, end):
                if b.strand=="-":
                    linestr1[x]="<";
                else:
                    linestr1[x]=">";

            # draw name
            endname = start + len(b.name)
            for x in range(start, endname):
                #print x, x-start
                linestr2[x]=b.name[x - start]

            # draw score
            if b.score!=0:
                scorestr = str(b.score)
                for x in range(start+1, start+1+len(scorestr)):
                    linestr1[x]=scorestr[x - start - 1]

        textLines.append( "".join(len(header) * [" "])+ linestr1.tostring() )
        textLines.append( "".join(len(header) * [" "])+ linestr2.tostring() )
    return "\n".join(textLines)

def prettyPrint_translated(seqs, allBeds, leftLimit=-1, rightLimit=-1):
    """ given already translated features and an alignment, format it in ascii format between leftLimit and rightLimit."""
    # put beds onto lines (y-movement calculation)
    seqlen = len(seqs[0].nucl)
    lines = putOntoLines(allBeds, seqlen)
    # convert beds and limits
    alnLimitLeft  = leftLimit
    alnLimitRight = rightLimit
    if leftLimit==-1: 
        alnLimitLeft = 0
    if rightLimit==-1:
        alnLimitRight = seqlen
    asciidata = paint(seqs, lines, alnLimitLeft, alnLimitRight) 
    return asciidata

        
def prettyPrintBlock(seqs, baseSeq, beds, leftLimit=-1, rightLimit=-1):
    #seqlen = len(baseSeq.nucl)
    # index beds by sequence name
    bedBySeq = {}
    for b in beds:
        bedBySeq.setdefault(b.chrom,[]).append(b)
    # index fasta seqs by seq name
    seqnameToIdx = {}
    for seq in seqs:
        if seq.id in seqnameToIdx:
            sys.stderr.write("error: fa id is not unique\n")
            sys.exit(1)
        seqnameToIdx[seq.id]=seq
        
    # put beds onto gapped sequence coordinates, gaps, translate window borders
    allBeds = []
    for seqName,beds in bedBySeq.iteritems():
        if seqName not in seqnameToIdx:
            sys.stderr.write("warning: cannot show features with seqname=%s, as there is no sequence for them\n" % seqName)
        else:
            seq = seqnameToIdx[seqName]
            (transBeds, alnLimitLeft, alnLimitRight) = ungappedToGapped(seq.nucl, beds, leftLimit, rightLimit)
            allBeds.extend(transBeds)

    (dummy, alnLimitLeft, alnLimitRight) = ungappedToGapped(baseSeq.nucl, [], leftLimit, rightLimit)
    if (alnLimitLeft, alnLimitRight == 0,0):
	alnLimitLeft, alnLimitRight = -1,-1
    asciidata = prettyPrint_translated(seqs, allBeds, alnLimitLeft, alnLimitRight)
    return asciidata


