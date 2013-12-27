from sys import *
import bed, util

# HELPER FUNCTIONS for all classes
def getClosestGenes(geneLeft, geneRight, distLeft, distRight, flankType, maxDist, stripComment, upstream):
    if stripComment:
        geneLeft = geneLeft.split("_")[0]
        geneRight = geneRight.split("_")[0]

    if geneLeft==geneRight:
        return [geneLeft], True

    genes = []

    if flankType==0: # only downstream genes
        if upstream==1:
            if distLeft < maxDist:
                genes.append(geneLeft)
        elif upstream==2:
            if distRight < maxDist:
                genes.append(geneRight)
        elif upstream==3:
            if distRight < maxDist:
                genes.append(geneRight)
            if distLeft < maxDist:
                genes.append(geneLeft)
        elif upstream==0:
            pass
        else:
            stderr.write("error parsing strands of flanking genes in wordfile\n")
            stderr.write("upstream parameter is"+str(upstream))
            exit(1)

    elif flankType==1: # the gene that is closer and < maxDist
        if distLeft < distRight and distLeft < maxDist:
            genes.append(geneLeft)
        elif distRight < maxDist:
            genes.append(geneRight)
    elif flankType==2: # both genes if they are closer than maxDist
        if distLeft < maxDist:
            genes.append(geneLeft)
        if distRight < maxDist:
            genes.append(geneRight)
    else:
        #stderr.write("internal program error: wrong flankType in wordParser.py\n")
        raise Exception('wrong flankType')

    return genes, False

class WordBlock:
    def __init__(self, chrom, start, end):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.words = {}
    
    def addWord(self, word, positions):
        """ positions is a list of tuples (pos, strand) """
        self.words[word] = positions

    def __repr__(self):
        return self.toString()

    def toString(self, number=None):
        lines = []
        lines.append("b\t%s\t%s\t%s" % (self.chrom, self.start, self.end))
        if number!=None:
            lines.append("n\t%d" % number)
        if "geneLeft" in self.__dict__:
            lines.append("f\t%s\t%d\t%s\t%d\t%d" % (self.geneLeft, self.distLeft, self.geneRight, self.distRight, self.upstream))
        lines.append("s\t%d" % (self.score))
        for w, posList in self.words.iteritems():
            lines.append("w\t%s\t%s," % (w, ",".join(["/".join(p) for p in  posList])))
        return "\n".join(lines)

    def getBeds(self, word):
        start = int(self.start)
        wordCount = len(self.words.get(word,[]))
        for p in self.words.get(word, []):
            if p == "":
                break
            #pos, strand = p.split("/")
            pos, strand = p
            pos = int(pos)
            yield(bed.Feature(self.chrom, start+pos, start+pos+len(word), word, wordCount, strand))
 
    def getGenes(self, flankType=1, maxDist=3000000000, stripComment=True):
        """ (genes, intronic) = return list of genes and wether the gene is intronic, if it's only one """
        """ flankType can be: 2=take both genes, 1=take closest gene, 0=only downstream genes (can lead to zero results) """
        if "geneLeft" in self.__dict__:
            geneLeft = self.geneLeft
            geneRight = self.geneRight
            distLeft = self.distLeft
            distRight = self.distRight
            return getClosestGenes(geneLeft, geneRight, distLeft, distRight, flankType, maxDist, stripComment, self.upstream)
        else:
            return [], False
        
    def getBedFeature(self):
        b = bed.Feature(seqid=self.chrom, start=self.start, end=self.end, score=self.score)
        return b

    def getCoverFeature(self):
        total = 0
        nameParts = []
        for w in self.words:
            count = sum([1 for b in self.getBeds(w)])
            total+=count
            nameParts.append(w+"="+str(count))
        b = bed.Feature(self.chrom, self.start, self.end, "_".join(nameParts), total, "+")
        return b
            
    def _checkSubsets(self, winSize, posList, wordCount):
        """
        we have len(wordCount) words
        wordCount contains the minimum number of words that we have to find
        posList contains a list of each word, with its start positions

        wordPos is out current pointer to posList, one per word

        currentWord is the index into wordPos that we are currently increasing

        algorithm:
            - check if the current window of wordCount positions in posList
              for all words is within winSize
            - if not: increase wordPos for currentWord 
        """
        #print "checkSubsets", winSize, posList, wordCount, "<p>"
        # determine maximum pos for each word
        posListMax = []
        for i in range(0, len(posList)):
            posListMax.append( len(posList[i])-wordCount[i])

        #print "posListMax", posListMax

        for wordPos in util.iterCounts(posListMax):
            # check subset represented by wordPos and return true if found
            pos = []
            for w in range(0, len(wordCount)):
                wps = posList[w][wordPos[w]:wordPos[w]+wordCount[w]]
                pos.extend([int(p) for p in wps])
            maxPos = max(pos)
            minPos = min(pos)
            diff = maxPos - minPos
            #print "pos", pos, "diff", diff, "winsize", winSize, "<p>"
            if diff < winSize:
                #print "FOUND in subroutine<P>"
                return True
        return False

    def _oneSubsetWithinWindow(self, words, b, winSize):
        """ static method: check if there is at least one subset of words in block that are closer together than winSize  """
        posList   = [] # a list of lists with positions for every word
        wordCount = [0] * len(words) # the number of words needed

        i = 0
        for w,count in words.iteritems():
            posList.append([p for p,s in b.words[w]])
            wordCount[i]=count
            i+=1

        #print "winsize", winSize, "<br>poslist",posList, "<br>wordcount", wordCount, "<br>wordpos", wordPos
        res = self._checkSubsets(winSize, posList, wordCount)
        #print "FOUND? %s<P>" % res
        return res

    def wordsFound(self,words, winSize=None):
        found=True
        for w in words:
            if w not in self.words:
                found=False
                #print "not found at all"
                break
            else:
                if not len(self.words[w]) >= words[w]:
                    found=False
                    #print "found, but not enough instances"
                    break
        if found:
            if not winSize or (winSize and self._oneSubsetWithinWindow(words, self, winSize)):
                #print "FOUND<P>"
                return True
        else:
            return False

class Words(list):
    def __init__(self):
        pass


def readMotifs(fh):
    """ retrieve list of motifs from words file, has to be called before readBlocks(), as the motifs are 
    at the start of the file. Returns [] if nothing found"""
    motifs = []
    firstLine = fh.readline()
    if not firstLine.startswith("##words"):
        stderr.write("wordParser.py/readWords: file %s is not a word-file\n" % fh.name)
        exit(1)

    l = fh.readline()
    while l!="":
        if l.startswith("b"):
            stderr.write("wordParser: no motifs found")
            exit(1)
        if l.startswith("m"):
            fs = l.strip().split()
            motif = fs[1].upper()
            desc = fs[2]
            motifs.append( (motif, desc) )
        if l=="\n":
            if len(motifs)!=0:
                return motifs
        l = fh.readline()
    return motifs

def readBlocks(fh):
    """ a GENERATOR function, can be used as an iterator"""
    block = None
    for l in fh:
        if l.startswith("b"):
            if block!=None:
                yield block
            b, chrom, start, end = l.strip().split()
            block = WordBlock(chrom, start,end)
            if fh.name=="<stdin>":
                block.pos = None
            else:
                block.pos = fh.tell()
        if l.startswith("n"):
            block.number=int(l.strip().split()[1])
        if l.startswith("f"):
            f, geneLeft, distLeft, geneRight, distRight, upstream = l.strip().split()
            block.geneLeft = geneLeft
            block.geneRight = geneRight
            block.distLeft = int(distLeft)
            block.distRight = int(distRight)
            block.upstream = int(upstream)
        if l.startswith("s"):
            s, score = l.strip().split()
            block.score = int(score)
        if l.startswith("d"):
            seq = l.strip().split()[5]
            if "seqs" not in block.__dict__:
                block.seqs = []
            block.seqs.append(seq)
        if l.startswith("w"):
            w, word, pos = l.strip().split()
            block.addWord(word, [p.split("/") for p in pos.split(",")[:-1]])
        if l.startswith("m"):
            pass
    if block!=None:
        yield block

def parseWords(strings):
    """ parse list of strings in the format GATTA:2 blalba:3 TTTNNTT into a dict """
    words = {}
    for a in strings:
        if a=="":
            continue
        if ":" in a:
            word, count = a.split(":")
        else:
            word = a
            count = 1
        words[word.upper()]=int(count)
    return words


# ==============BLOCK PARSER================
# functions to parse index files for blocks
def parseIndex(fh):
    blockToGene = {}
    motifToBlockOccs = {}

    for l in fh:
        # example
        # bf	0	ci0100135030_exon6	18	ci0100135008_exon	2782	1
        if l.startswith("bf"):
            bf, no, geneLeft, distLeft, geneRight, distRight, upstream = l.strip().split()
            distLeft = int(distLeft)
            distRight = int(distRight)
            maxDist = 300000000
            flankType = 1 # only closest gene
            genes, intronic = getClosestGenes(geneLeft, geneRight, distLeft, distRight, flankType, maxDist, True)
            blockToGene[int(no)]=genes[0]
        else:
            #example CGCNNCGC	3073/1,5122/1,10253/1,1806/1,2449/1,8082/1,8597/1,12567/1,4633/1,12443/1,286/1,48/1,10914/1,6947/1
            motif, occs = l.strip().split("\t")
            assert(not motif in motifToBlockOccs)
            occs = occs.split(",")
            motifOccs = {}
            for occPart in occs:
                if occPart=="":
                    continue
                block, occ = occPart.split("/")
                block=int(block)
                occ = int(occ)
                motifOccs[block]=occ
            motifToBlockOccs[motif]=motifOccs
    return motifToBlockOccs.keys(), blockToGene, motifToBlockOccs
