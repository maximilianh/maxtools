import sys, os, textwrap

class FastaSeq:
    """ parse UCSC-style Fasta files with special fasta headers 
     Example: >hg17 range=chr1:1-50 
     """
    def __init__(self, id, seq):
        self.id = id.split()[0]
        self.fullid = id
        self.completeid = id # for compat.
        self.nucl = seq 
        try:
            range = id.split(" ")[1].split("=")[1]
            chromcoord = range.split(":")
            self.chrom = chromcoord[0]
            coords = chromcoord[1].strip().split("-")
            self.start = int(coords[0])
            self.end = int(coords[1])

            self.strand=None
            strandTok = id.split(" ")[2]
            strandToks = strandTok.split("=")
            if strandToks[0].lower()=="revcomp":
                if strandToks[1].lower()=="false":
                    self.strand="+"
                if strandToks[1].lower()=="true":
                    self.strand="-"

        except IndexError:
            #sys.stderr.write("Just for info: range=x:y-z not found in fasta header of sequence with id '"+id+"' . \n")
            self.chrom = self.id
            self.start = 0
            self.end   = len(self)
    def __len__(self):
        return len(self.nucl)
    def __repr__(self):
        return self.toFasta()

    def toFasta(self, ranges=True):
        str = ""
        if ranges:
            str += ">"+self.id+" range=%s:%d-%d\n" % (self.chrom, self.start, self.end)
        else:
            str += ">"+self.fullid+"\n"
        str += self.nucl
        return str

    def getId(self):
        return self.id

    def getChrom(self):
        return self.chrom

    def getID(self):
        return self.id.split()[0]

    def getStart(self):
        return self.start

    def getEnd(self):
        return self.end

    def getNucl(self):
        return self.nucl
    def getUngapNucl(self):
        return self.nucl.replace("-","")

    def getIdField(self, i):
        return self.fullid.split("|")[i]

    def getEnsemblGene(self):
        return self.getIdField(3)

    def getEnsemblTranscript(self):
        return self.getIdField(4)

def readFastaLines(lines):
  seqs = {}
  seqLines = []
  id = ""
  for line in lines:
      if line.startswith(">"):
            if seqLines!=[]: # on first >, seq is empty
                  faseq = FastaSeq(id, "".join(seqLines))
                  seqs[id.split()[0]]=faseq
                  seqLines=[]
            id=line.strip(">").strip()
      else:
          #seqLines.append(line.replace(" ","").strip())
          seqLines.append(line.strip())
  seqs[id.split()[0]]=FastaSeq(id, "".join(seqLines)) # add last seq
  return seqs

def readFastaFh(f):
  """ parses fasta from filehandle """
  lines=f.readlines()
  return readFastaLines(lines)

def readFasta(filename):
  ## will read from lines
  ## will return a dict of seqid -> FastaSeq()
  if filename!="stdin":
      f = open(filename,"r")
  else:
      f = sys.stdin
  return readFastaFh(f)

def readFastaAsList(filename):
  ## will read from lines
  ## will return a  list of seqs
  if filename!="stdin":
      f = open(filename,"r")
  else:
      f = sys.stdin

  lines=f.readlines()
  seqs = []
  seq = ""
  id = ""
  for line in lines:
          if line.startswith("#") or line.startswith("="):
              continue
          if line.startswith(">"):
            if seq!="": # on first >, seq is empty
                  faseq = FastaSeq(id, seq)
                  seqs.append(faseq)
                  seq = ""
            id=line.strip(">").strip()
          else:
            seq+=line.replace(" ","").strip()
  if seq=="":
      return None
  else:
      seqs.append(FastaSeq(id, seq) )
      return seqs

def parseFastaUntil(f, max, lastline):
  """ reads seqs into memory from filehandle f, up to max, inf. if max=-1, optional oldid is the header line and will be parsed """
  """ last id line has to be read in by caller and supplied, and provided at every call 
      returns tuple: """
  seqs = []
  lines = []
  id = lastline.strip(">").strip()
  count = 1
  stopped = False
  lastline = ""

  for line in f:
          if line.startswith("#"):
              continue
          if line.startswith(">"):
                count += 1
                if count > max:
                    lastline = line
                    stopped=True
                    break
                if len(lines)!=0: # on first >, seq is empty
                      faseq = FastaSeq(id, "".join(lines))
                      seqs.append(faseq)
                      id=line.strip(">").strip()
                      lines=[]
          else:
                lines.append(line.replace(" ","").strip())
  if len(lines)!=0:
      faseq = FastaSeq(id, "".join(lines))
      seqs.append(faseq)
  return lastline, seqs


# ======================
def readMaf(filename):
    """ parses a maf file and will return a list (blocks) of 
        a dict (block) of fasta-sequences """

    f = open(filename, "r")
    block = {}
    blocks = []
    blockno = 1

    for l in f:
        if l.startswith("#"):
            continue
        if l.startswith("a "):
            if block!={}:
                blocks.append(block)
                block = {}
                blockno += 1
            score = float(l.split("=")[1].strip())
        if l.startswith("s "):
            parts = l.split()
            if "." in parts[1]:
                org, chrom = parts[1].split(".")
            else:
                org = chrom = parts[1]
            start = int(parts[2])
            end = start+int(parts[3])
            strand = parts[4]
            nucl = parts[6].strip()
            seq = FastaSeq("%s range=%s:%u-%u origin=%s block=%u score=%3.f strand=%s" \
                    % (org, chrom, start, end, filename, blockno, score, strand), nucl)
#            if block.get(org)==None:
#                block[org]=[]
            block[org]=seq
    blocks.append(block)
    return blocks

def mafToFasta(block):
    """ converts one single maf-block to a list of fasta sequences"""
    faseqs = []
    for org, fasta in block.iteritems():
        faseqs.append(fasta)
    return faseqs

def mafBaseOrg(fname):
    """ return base species in maf file """
    for l in open(fname, "r"):
        if l.startswith("s"):
            fs = l.split()
            org = fs[1].split(".")[0]
            break
    return org

# ========================
def readMafNext(file):
    """ parses a maf file from filehandle and will return the next block of aligned fasta-sequences as a dict orgname -> FastaSequence. Returns None if EOF. """

    block = {}
    blocks = []
    blockno = 1

    for l in file:
        if l.startswith("#"):
            continue
        if l.startswith("a "):
            if block!={}:
                return block
        if l.startswith("s "):
            parts = l.split()
            org, chrom = parts[1].split(".")
            start = int(parts[2])
            end = start+int(parts[3])
            strand = parts[4]
            nucl = parts[6].strip()
            seq = FastaSeq("%s range=%s:%u-%u origin=%s block=%u strand=%s" \
                    % (org, chrom, start, end, file.name, blockno, strand), nucl)
            block[org]=seq
    return block

# ===========================
def readFastaa(fname):
    """ ensembl format: alignments in fasta format, separated by # """
    sys.stderr.write("reading %s\n" % fname)
    fh = open(fname, "r")
    allseqs = []
    lastline=fh.readline()
    while True:
        lastline, seqs = parseFastaUntil(fh, 2, lastline)
        if len(seqs)==0:
            break
        allseqs.append(seqs)
    return allseqs

# ==========================
def fastaIdReformat(outFormat, inFileName, of, inplace):
    """ change the id lines of a fasta file, inplace will modify the file itself, of is outfile handle """
    if inFileName=="stdin":
        fh = stdin
    else:
        fh = open(inFileName)
        if inplace:
            tmpout = inFileName+".tmp"
            of = open(tmpout, "w")

    for l in fh:
        if not l.startswith(">"):
            of.write(l)
        else:
            id = l.strip().strip(">")

            if outFormat=="ensTrans":
                # >ENST00000400685 cdna:novel supercontig::NT_113903:9607:12778:1 gene:ENSG00000215618
                parts = id.split(" ")
                gene = parts[-1].strip("gene:")
                trans = parts[0]
                # remove Ensembl-Genes prefix
                trans = trans.replace("EG:", "")
                id = trans+"|"+gene
            if outFormat=="strip":
                id = id.split()[0]
                id = id.split("#")[0]
            else:
                stderr.write("error: outFormat %s not recognized\n")
                exit(1)

            of.write(">"+id+"\n")

    if inplace:
        os.rename(tmpout, inFileName)
