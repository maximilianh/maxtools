# tools to work with bed files

from sys import *
from optparse import OptionParser
from bed import Feature,Features,parseBedFilename
import os
import tempfile

QUIET=False


# ====== BedFindNeighbor ======

def debug(msg):
    if not QUIET:
        stderr.write(msg+"\n")

gapFeature = Feature("noChrom",0,0,"Gap",999,".")
n = 0 

def isExon(bed):
    return bed.name.endswith("#")

def searchNextExon(bedList,i):
    if i>=len(bedList)-1:
        return gapFeature
    #print "search, i=", i
    #print "search:", bedList[i].name
    while not isExon(bedList[i]):
        #print "check:", bedList[i].name
        i+=1
        if i==len(bedList)-1:
            return gapFeature
    #print "return:", bedList[i].name
    return bedList[i]

def createFeatures(bed, exonLeft, exonRight, noStar, addDist, onlyDownstream, onlyClosest):
    """ for bedFindNeighbors """
    def createFeature(bed, name, addDist, dist=None):
        if addDist:
            if not dist: # add two distances
                distLeft = bed.start - exonLeft.end
                distRight = exonRight.start - bed.end
                name="%s|%d|%d" % (name, distLeft, distRight)
            else: # add only one distance 
                name = "%s|%d" % (name, (-1)*dist) # biologists prefer negative values in the upstream region

        #cols =  [bed.chrom, str(bed.start), str(bed.end), bed.name+"|"+name, str(bed.score), bed.strand]
        #print "\t".join(cols)
        newName= bed.name+"|"+name
        b = Feature(bed.chrom, bed.start, bed.end, newName, bed.score, bed.strand)
        global n
        n+=1
        return b

    beds = Features()
    if exonLeft.chrom!=bed.chrom:
        exonLeft=gapFeature
    if exonRight.chrom!=bed.chrom:
        exonRight=gapFeature
    name1 = exonLeft.name.strip("#")
    name2 = exonRight.name.strip("#")
    # -> always two features, no stars
    if noStar: 
        beds.append(createFeature(bed, name1+"|"+name2, addDist))
    # -> generate 0-2 features (only if gene is downstream)
    elif onlyDownstream: 
        if exonLeft.strand=="-":
            beds.append(createFeature(bed, name1, addDist, bed.start - exonLeft.end))
        if exonRight.strand=="+":
            beds.append(createFeature(bed, name2, addDist, exonRight.start - bed.end))
    # -> always generate one feature
    elif onlyClosest: 
        if exonLeft.strand=="+":
            leftStart=exonLeft.start
        else:
            leftStart =exonLeft.end
        if exonRight.strand=="+":
            rightStart=exonRight.start
        else:
            rightStart =exonRight.end
        distLeft = abs(bed.start - leftStart)
        distRight = abs(rightStart - bed.end)
        if distLeft < distRight:
            closeName = name1
        else:
            closeName = name2
        beds.append(createFeature(bed, closeName, addDist))

    # default, always two features, with stars
    else: 
        if exonLeft.strand=="-":
            name1 = "*"+name1
        if exonRight.strand=="+":
            name2 = "*"+name2
        beds.append(createFeature(bed, name1+"|"+name2, addDist))
    return beds

def bedFindNeighbors(regions,genes,noStar,onlyDownstream,onlyClosest,addDist,keepOverlaps):
    tmpRegFh, tempRegFname = None,None
    if type(regions)==type([]) or type(regions)==type(Features()):
        tmpRegFh, tempRegFname = tempfile.mkstemp()
        regions.writeToFile(tempRegFname)
        regions = tempRegFname

    tmpfname = tempfile.mktemp()
    debug("concatting and sorting bed files %s and %s into %s...\n" % (genes, regions, tmpfname))
    # add marker to bed features and sort
    cmd = 'cat %s | gawk \'{OFS="\t"; $4=$4"#"; print }\' | cat - %s | bedSort stdin stdout > %s' % (genes, regions, tmpfname)
    debug("exec: %s\n" % cmd)

    ret = os.system(cmd )
    if ret != 0:
        stderr.write("Error executing command: %s" % cmd)
        exit(1)

    debug("Reading temp bed file...\n")
    beds = parseBedFilename(tmpfname, False, QUIET)

    lastExon = gapFeature
    nextExon = gapFeature

    debug("Annotating features and printing...\n")
    newbeds = Features()
    # iterate over the mixed bed features
    for i in xrange(0, len(beds)):
        b = beds[i]
        if isExon(b):
            lastExon=b
        else:
            #print "lastExon, nextExon", lastExon.name, nextExon.name
            if nextExon==lastExon or nextExon.start < lastExon.start or lastExon.chrom!=nextExon.chrom:
                nextExon = searchNextExon(beds, i+1)
            if b.overlaps(lastExon):
                debug("info: region %s is overlapping (last) exon %s\n" % ( b, lastExon))
                if not keepOverlaps:
                    continue
            if b.overlaps(nextExon):
                debug("info: region %s is overlapping (next) exon %s\n" % ( b, nextExon))
                if not keepOverlaps:
                    continue
            newbs =  createFeatures(b, lastExon, nextExon, noStar, addDist, onlyDownstream, onlyClosest)
            newbeds.extend(newbs)

    os.remove(tmpfname)
    if tempRegFname!=None:
        os.remove(tempRegFname)
    return newbeds
    #debug("Obtained %d features.\n" % len(newbeds))

