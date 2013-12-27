class SeqWindow:
    """ a window on a gapped sequence that can be shifted to the right. It has a flexible size, but
        will always contain w characters. We use UCSC style coordinates.
    """
    def __init__(self, seq, windowSize):
	self.seq=seq
	self.w = windowSize
	(start,end) = self.initWindow(windowSize)
	self.start = start 
	self.end = end
        self.endReached = False

    def nucl(self):
        """ the nucleotide sequence of the window """
        self.seq.nucl[self.start:self.end].replace("-","")

    def searchNonGap(self, pos):
	""" returns the first non-gap position after pos. Return 0 if exceeded length. """
	if pos+1>=len(self.seq.nucl):
	    return 0
	pos+=1
	while self.seq.nucl[pos]=="-":
	    if pos+1>=len(self.seq.nucl):
		return 0
	    pos+=1
	return pos

    def stepUntilXChars(self, pos, w):
	""" returns the next position p>pos such that there are at least w non-gap chars between pos and p. Return 0 if sequence is too short. """
	#gaps = seq[pos,pos+x+1].count("-")
	endPos = pos
	chars=0
	while chars<w:
	    endPos = self.searchNonGap(endPos)
            if endPos==0:
                return 0
	    chars+=1
	return endPos

    def initWindow(self, w):
	""" will return (start,end) coordinates such that there are at least w non-gaps, starting from 0. Will return (0,0) if seq is too short.  """ 
	start = 0
        #print type(self.seq)
        #print self.seq
	if self.seq.nucl[0]=="-":
	    start = self.searchNonGap(0)
	end = self.stepUntilXChars(start, w)
        #print " init ","seq="+self.seq.id, start
        ## TODO THIS WILL STOP IF len(seq) == len(window) !!!
        if end==0:
            raise StopIteration
        else:
            return (start,end)

    def shiftRight(self):
	""" shifts right one position """	
	newEnd = self.searchNonGap(self.end-1)
	if newEnd==0:
	    self.endReached=True
	    return
	else:
	    self.end=newEnd+1
	    self.start = self.searchNonGap(self.start)

