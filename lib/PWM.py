import re
from math import log
import copy
import SubSeq
import bed
import sys

class pwms(dict):
    """ this is returned by transfac parsers """
    " A list of pwms "
    def __repr__(self):
        lines = []
        for p in self.values():
            #lines.append(p.ac+"\t"+p.id+"\t"+str(p.getRealIC())+"\t"+str(len(p.siteseqs)))
            lines.append(p.ac+"\t"+p.id+"\t"+str(len(p.siteseqs)))
        return "\n".join(lines)

class pwm:
   """ matrix is a list (for every pos) of dictionaries (keys: a,c,t,g)
    matrix has additional attributes. you can scan a sequence against a pwm."""
   def __init__(self):
        self.ac = ""
        self.id = ""
        self.description = ""
        self.factors = []
        self.factoracs = []
        self.siteacs = []
        self.siteseqs = []
        self.sitestrands = []
        self.sitestarts = []
        self.sitelens = []
        self.weights = []
        self.float = False # note if the original matrix consistet of floats
        # this is only to speed up calculations 
        self.minscore=0
        self.maxscore =0
        self.IC = 0
        self.cutoff = 0
        self.mostConserved = []
        self.coreMinscore = 0
        self.coreMaxscore = 0
        self.coreLookup = {}
   def __len__(self):
        return len(self.weights)
   def _regularize(self):
        """ convert counts to probabilities """
        newmatrix = []
        counts = self.weights
        for v in counts:
          total = sum(v.values())
          # count -> freq
          freq={}
          for nucl,count in v.iteritems():
                freq[nucl]=float(count/total)
          newmatrix.append(freq)
        self.weights = newmatrix
   def _pseudocounts(self):
        """ add pseudocounts to matrix """
        newmatrix = []
        for w in self.weights:
          freq={}
          for nucl,prob in w.iteritems():
                freq[nucl]=prob+0.001
          newmatrix.append(freq)
        self.weights=newmatrix
   def _mostConserved(self, x):
       """ will create a list of the x most conserved nucleotides """
       #get x highest positions
       l = []
       for i in range(0, len(self.weights)):
           (pos, IC) = (i, self.IC[i])
           l.append( (pos,IC) )
       l.sort(key=lambda (pos, IC): IC, reverse=True) 
       l = l[0:x] 
       pos = [x for (x,y) in l]
       pos.sort()
       pos = set(pos)
       self.mostConserved = pos
       
   def prepareScoring(self):
       """ this will regularize and change to matrix to make scanning faster.
       You have to call this before you calling any scoring functions """
       self._regularize()
       self._pseudocounts()
       self._calcMatchIC()
       self._mostConserved(5)
       self.minscore, self.maxscore = self._match_calcminmax(range(0, len(self.weights)), self.IC)
       self.coreMinscore, self.coreMaxscore = self._match_calcminmax(self.mostConserved, self.IC)
       self._calcCoreScoreLookup()
                    
   def _calcCoreScore(self,nucl):
        """ calc corescore for given core (=nucl= string of length 5) """
        current=0
        corePos=0
        # iterate over mostConservPosition, for each: get Probab, extract p for nucl, extract IC for nucl
        for pos in self.mostConserved:
            probs = self.weights[pos]
            nuc = nucl[corePos].lower()
            corePos+=1
            IC = self.IC[pos]
            current += probs[nuc]*IC
        #print self.coreMaxscore, self.coreMinscore
        corescore = (current - self.coreMinscore) / (self.coreMaxscore - self.coreMinscore)
        #print nucl, corescore
        return corescore

   def _calcCoreScoreLookup(self):
        """ this is building a dict of a dict of dict of a dict of 
            every possible core 4mer -> core score to make this part 
            of the program faster and avoid calc. of logs """
        nucl = ['a', 'c', 'g', 't']
        table = {}
        for a1 in nucl:
            table[a1] = {}
            for a2 in nucl:
                table[a1][a2] = {}
                for a3 in nucl:
                    table [a1][a2][a3] = {}
                    for a4 in nucl:
                        table[a1][a2][a3][a4] = {}
                        for a5 in nucl:
                            table [a1][a2][a3][a4][a5]  = self._calcCoreScore([a1,a2,a3,a4,a5])
        self.coreLookup = table

   def __repr__(self):
       """ output in transfac format """
       print "AC  ",self.ac
       print "ID  ",self.id
       print "DE  ",self.description
       i = 1
       lines = []
       lines.append( "P0     A       C      G     T")
       for w in self.weights:
           #lines.append("%02d   %.3f  %.3f  %.3f  %.3f" % (i, w['a'],w['c'],w['g'],w['t']))
           lines.append("%02d   %u  %u  %u  %u" % (i, w['a'],w['c'],w['g'],w['t']))
           i += 1
       lines.append("XX")
       lines.append("CC  Cutoff is set to: "+str(self.cutoff))
       lines.append("CC  Five most conserved positions: "+str([x+1 for x in self.mostConserved]))
       lines.append( "//")
       return "\n".join(lines)

   def _calcMatchIC(self):
        # you need to calc pseudo counts before calling this function!
        # calculates the Ci-vector as explained
        # in Kel 03, a value in 0..100
        # for every pos in the matrix that describes 
        # the conservation of a position
        # and is used to scale the pos-scores.
        C = []
        # print mat
        for probs in self.weights:
                ics = [(p * log(4*p)) for p in probs.values()]        
                C.append(sum(ics))
        self.IC=C

   def log2(self, x):
       """ calc logarithmus dualis """
       return log(x) / log (2)

   def getRealIC(self):
        # calc 2 - sum(p * ld (p)) and sum over all positions
        C = []
        # print mat
        for probs in self.weights:
                ics = [(p * self.log2(p)) for p in probs.values()]
                C.append(2 + sum(ics))
        return sum(C)

   def _match_calcminmax(self, positions,cons):
        ### calc minimum/maximum possible score for a matrix
        ### see match paper
        maxscore = float(0)
        minscore = float(0)
        for pos in positions:
                probs = self.weights[pos]
                IC = self.IC[pos]
                maxscore+=max(probs.values())*IC
                minscore+=min(probs.values())*IC
        return (minscore, maxscore)
                    
   def match(self, subseq, corecutoff, cutoff):
       """ annotates a subseq with tfname, score and ref to pwm if it matches """
       (corescore, matscore) = self.score(subseq.getUngapNucl(), corecutoff) 
       if not (corescore >= corecutoff and matscore >= cutoff):
               return None
       else:
           subseq.name = self.id
           subseq.score = matscore
           subseq.pwm = self
           return subseq

   def score(self, nucl, corecutoff):
        """ score matrix against sequence with Biobase Match (Kel et al 2003) method. return 2-tuple of
        corescore, matrixscore or corescore,0 if corescore did not exceed corecutoff.
        Is using the coreLookup table for fast log-calculation.
        """
        if not len(self.weights)==len(nucl):
            raise Exception, "matrix length has to be the same as sequence length. len(matrix)="+str(len(self.weights))+",len(seq)="+str(len(nucl))
        # calc socre on just core
        current = 0
        core = []
        for pos in self.mostConserved:
            core.append(nucl[pos].lower())
        #        if "n" in core:
        #            sys.stderr.write("ERROR! FOUND AN N!")
        #            print nucl
        corescore = self.coreLookup [core[0]] [core[1]] [core[2]] [core[3]] [core[4]]
        #print corescore

        if corescore<corecutoff:
            return (corescore, 0)
        # calc on whole matrix
        current = 0
        for pos in xrange(0, len(nucl)):
                nuc = nucl[pos].lower()
                probs = self.weights[pos]
                IC = self.IC[pos]
                current += probs[nuc]*IC
        score = (current - self.minscore) / (self.maxscore - self.minscore)
        return (corescore, score)

   def isInCore(self, pos):
       """ checks if a give nucleotide is part of the x most conserved positions """
       pass

def readMatrices(filename):
    """ will read matrices from transfac-style file returns a dict of matrixid -> pwm 
    example matrix: 
    ID TEST
    P0 A C T G
    01 1 2 0 2
    02 1 2 3 0 """
    if filename!="stdin":
        lines = open(filename, "r").readlines()
    else:
        lines = sys.stdin.readlines()

    matrices = pwms()
    start = False
    matrix = pwm()
 
    # read data into matrices
    for line in lines:
      line = line.strip()
      if line.startswith("ID"):
              id = line.split()[1].strip()
              matrix = pwm()
              matrix.id = id
              matrix.ac = lastacc
              starts = False
              continue
      if line.startswith("AC"):
              ac = line[2:].strip()
              lastacc = ac
      if line.startswith("P0"):
            start = True
            continue
      if line.startswith("//") and matrix.weights != []:
              if matrix.id=="":
                      sys.stderr.write("error: no id for matrix")
              matrices[matrix.id] = matrix
      if line.startswith("BF"):
              data = line[2:].strip().split(";")
              factorac = data[0]
              factor = data[1]
              matrix.factoracs.append(factorac)
              matrix.factors.append(factor)
      if (line.startswith("DE")):
              matrix.description = line[2:].strip()
      if (line.startswith("BS")):
              rest = line[2:].strip()
              parts = rest.split(";")
              siteseq = parts[0].strip()
              siteac = parts[1].strip()
              start = parts[2].strip()
              len = parts[3].strip()
              strand = parts[5].strip()
              if strand=="p.":
                  strand = "+"
              else:
                  strand = "-"
              matrix.siteacs.append(siteac) 
              matrix.siteseqs.append(siteseq)
              matrix.sitestrands.append(strand)
              matrix.sitestarts.append(start)
              matrix.sitelens.append(len)
      # OK, found start of matrix
      if start and re.compile("^[0-9][0-9] ").match(line):
              counts = line.split()[1:]
              counts = counts[:4]

              for c in counts:
                  if "." in c:
                      matrix.float = True
                  else:
                      matrix.float = False
              counts = [float(x) for x in counts ]
              probs = {'a': counts[0], 'c':counts[1], 'g':counts[2], 't':counts[3]}

              counts = [float(x) for x in counts ]
              matrix.weights.append(probs)
    return matrices
