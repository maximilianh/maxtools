import re
from math import log
import copy
import SubSeq
import bed

class pwms(dict):
    """ this is returned by transfac parsers """
    " A list of pwms "
    def __repr__(self):
	lines = []
	for p in self.values():
	    lines.append(p.ac+"\t"+p.id)
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
       """ will make a list of the x most conserved nucleotides """
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

   def __repr__(self):
       """ output in transfac format """
       print "AC  ",self.ac
       print "ID  ",self.id
       print "DE  ",self.description
       i = 1
       lines = []
       lines.append( "P0     A       C      G     T")
       for w in self.weights:
	   lines.append("%02d   %.3f  %.3f  %.3f  %.3f" % (i, w['a'],w['c'],w['g'],w['t']))
	   i += 1
       lines.append("CC  Cutoff is set to: "+str(self.cutoff))
       lines.append("CC  Five most conserved positions: "+str([x+1 for x in self.mostConserved]))
       lines.append( "//")
       return "\n".join(lines)
   def _calcMatchIC(self):
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
       (corescore, matscore) = self.score(subseq.getUngapNucl(), corecutoff) 
       if not (corescore >= corecutoff and matscore >= cutoff):
               return None
       else:
#           if self.id=="V$FOXD3_01":
#               print subseq,
#               print " ",corescore, matscore
           subseq.name = self.id
           subseq.score = matscore
           subseq.pwm = self
           return subseq

   def score(self, nucl, corecutoff):
        """ score matrix against sequence with Biobase Match (Kel et al 2003) method. return 2-tuple of
        corescore, matrixscore or corescore,0 if corescore did not exceed corecutoff.
        """
	if not len(self.weights)==len(nucl):
	    raise Exception, "matrix length has to be the same as sequence length. len(matrix)="+str(len(self.weights))+",len(seq)="+str(len(nucl))
        # just core
        current = 0
        for pos in self.mostConserved:
		nuc = nucl[pos].lower()
		probs = self.weights[pos]
		IC = self.IC[pos]
		current += probs[nuc]*IC
        #print self.coreMaxscore, self.coreMinscore
	corescore = (current - self.coreMinscore) / (self.coreMaxscore - self.coreMinscore)
        if corescore<corecutoff:
            return (corescore, 0)
        # whole matrix
	current = 0
	for pos in range(0, len(nucl)):
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
    lines = open(filename, "r").readlines()
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
	      starts = False
	      continue
      if line.startswith("AC"):
	      ac = line[2:].strip()
	      matrix.ac = ac
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
	      siteac = line[2:].strip().split(";")[1].strip()
              matrix.siteacs.append(siteac) 
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
