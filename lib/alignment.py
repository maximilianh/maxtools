import Fasta


class Alignment(list):
    """ a list of FastaSeq """
    def __init__(self, seqList=[]):
        for s in seqList:
            self.append(s)

    def __repr__(self):
        lines = []
        lines.append("Alignment of "+str(len(self))+" sequences")
        for s in self:
            lines.append(s.id+":"+s.nucl[0:20]+"..."+s.nucl[-20:])
        return "\n".join(lines)

    def writeToFasta(self, filename, gaps, ranges=True):
        blocks = self
        file = open(filename, "w")
        if gaps:
            file.write(self.asFasta(ranges))
        else:
            file.write(self.asFastaNoGaps())
        file.close()

    def asFasta(self, ranges=True):
        seqs = []
        for seq in self:
            seqs.append(seq.toFasta(ranges))
        return "\n".join(seqs)+"\n"

    def asFastaNoGaps(self, idPrefix=""):
        seqs = []
        for seq in self:
            seqs.append(">"+str(idPrefix)+str(seq.id)+"\n"+str(seq.nucl.replace("-","")))
        return "\n".join(seqs)+"\n"


