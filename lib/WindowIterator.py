import Fasta
from SeqWindow import *
from SubSeq import *

def subseq_generator(seq, winsize):
    """ returns gapless subsequences of fixed size from a gapped sequence """
    win = SeqWindow(seq, winsize)
    while not win.endReached:
        yield SubSeq(seq, win.start, win.end, "+" )
        yield SubSeq(seq, win.start, win.end, "-" )
        win.shiftRight()



