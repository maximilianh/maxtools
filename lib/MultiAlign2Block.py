from SeqWindow import *

class Block(list):
    """ a list of subsequences """
    def __repr__(self):
        str = ""
        for i in self:
            str += repr(i) 
        return str

    def getBlockPos(self, seqid):
        """ return absolute start of block on sequence with given seqid """ 
        ss = self.searchSeq(seqid)
        return ss.getAbsUngapStart()
    
    def searchSeq(self, seqid):
        """ return sequence with given seqid """
        found = None
        for ss in self:
            if ss.seq.getID()==seqid:
                found = ss 
                break
        if found==None:
            raise Exception, "A sequence with seqid '"+seqid+"' does not exist in alignment"
        return found

    def getBedLine(self, org):
        """ return a string in bed-format with ungapped coords of subseq from org """
        found = self.searchSeq(org)
        return found.seq.getChrom()+" "+str(found.getAbsUngapStart())+" "+str(found.getAbsUngapEnd())

class WindowIterator:
    """an iterator the will return subsequences of fixed sizes on longer sequences """
    def __init__(self, seq, winsize):
        self.seq = seq
        self.pos = 0
        self.winsize=winsize
        self.lastpos = len(self.seq.nucl)-self.winsize+1
    def __iter__(self):
        return self
    def next(self):
        if self.pos  == self.lastpos:
            raise StopIteration
        else:
            ss = Subseq(self.seq,self.pos,self.pos+self.winsize)
            self.pos +=1
            return ss

# --------


