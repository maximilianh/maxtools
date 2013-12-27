import Fasta

def revcomp(seq):
    table = { "a" : "t", "A" : "T", "t" : "a", "T" : "A", "c" : "g", "C":"G", "g":"c", "G":"C", "-":"-", "N":"N", "n":"n" }
    newseq = []
    for nucl in reversed(seq):
        newseq.append(table[nucl])
    return newseq

class SubSeq:
    """ a short, gapped, stranded feature located on another sequence. Can have a name. """
    def __init__(self, seq, start, end, strand):
        self.seq = seq
        self.start = start
        self.end = end
        self.gaplesslen = len(self.seq.getUngapNucl())
        self.strand = strand
        self.name = None
        self.score = None
    def __len__(self):
        return self.end-self.start
    def getAbsUngapStart(self):
        """ returns start pos on parent sequence, ungapped """
        return self.start-self.seq.nucl.count("-",0,self.start)+self.seq.getStart()
    def getAbsUngapEnd(self):
        return self.end-self.seq.nucl.count("-",0,self.end)+self.seq.getStart()
    def getAbsStart(self):
        """ returns start pos on parent sequence """
        return self.start+self.seq.getStart()
    def getAbsEnd(self):
        """ returns end pos on parent sequence """
        return self.end+self.seq.getStart()
    def getUngapNucl(self):
        return [x for x in self.getNucl() if x !="-"]
    def getNucl(self):
        seq = self.seq.getNucl()[self.start:self.end]
        if self.strand=="-":
            return revcomp(seq)
        else:
            return seq
    def getChrom(self):
        return self.seq.getChrom()
    def __repr__(self):
        return self.seq.getID()+" "+self.seq.getChrom()+" "+str(self.getAbsUngapStart())+" "+str(self.getAbsUngapEnd())+" "+self.strand+" "+"".join(self.getUngapNucl())

    def getBedLine(self, genomeCoords):
        if genomeCoords:
            return "%s\t%u\t%u\t%s\t%u\t%s" % (self.seq.getChrom(), self.getAbsUngapStart(), self.getAbsUngapEnd(), self.name, self.score*1000, self.strand)
        else:
            return "%s\t%s\t%s\t%s\t%u\t%s" % (self.seq.id, self.start, self.end, self.name, self.score*1000, self.strand)

    def getUngapLength(self):
        return self.gaplesslen

    def conservedPos(self):
        """ return positions with a high IC of matching pwm in absolute coords """
        start = self.getAbsUngapStart()
        pos = self.pwm.mostConserved
        pos = [x+start for x in pos]
        pos = set(pos)
        return pos

