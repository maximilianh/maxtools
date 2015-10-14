# search sequences for IUPAC motifs

import sys
import re
# jit compiler for python, remove these lines if you get an error here, it is no
try:
    import psyco
except (ImportError, ):
    pass
else:
    psyco.full()

# sorted alphabetically
nuclToIupac = {
    "CT"   : "Y",
    "AG"   : "R", 
    "AC"   : "M",
    "GT"   : "K",
    "CG"   : "S",
    "AT"   : "W",
    "ACT"  : "H",
    "CGT"  : "B",
    "ACG"  : "V",
    "AGT"  : "D",
    "ACGT" : "N",
    }


iupacTable = { "Y" : "TCY", "R" : "GAR", "M" : "ACM", "K" : "GTK", "S" : "GCS", "W" : "ATW", "H" : "ACTHYKW", "B" : "GTCBKYS", "V" : "GCAVSR", "D" : "GATDRWK", "N" : "ACTGNYRMKWSHBVD",
            "y" : "tcy", "r" : "gar", "m" : "acm", "k" : "gtk", "s" : "gcs", "w" : "atw", "h" : "acthykw", "b" : "gtcbkys", "v" : "gcavsr", "d" : "gatdrwk", "n" : "actgnyrmkwshbvd"}

def resolveIupac(seq):
    """ convert string with iupac characters to a regular expression string """
    newseq = []
    error=False
    for nucl in seq:
       if nucl in iupacTable:
           newseq += "[%s]" % iupacTable[nucl]
       else:
           newseq += nucl
           error=True
    #if error:
        #sys.stderr.write("warning: %s does not look like a iupac string\n" % seq)
    newstr = "".join(newseq)
    #newstr = newstr.replace("N", "[ACTGN]")
    return newstr

def findRegex(regEx, string, posList, strand):
        iter = regEx.finditer(string)
        for m in iter:
            #print m.start()
            posList.append((m.start(),strand))

def revComp(seq):
    table = { "a":"t", "A":"T", "t" :"a", "T":"A", "c":"g", "C":"G", "g":"c", "G":"C", "N":"N", "n":"n", 
            "Y":"R", "R" : "Y", "M" : "K", "K" : "M", "W":"W", "S":"S",
            "H":"D", "B":"V", "V":"B", "D":"H" }
    newseq = []
    # for python 2.3
    seq = list(seq)  # string to list
    seq.reverse()     
    # end python 2.3
    for nucl in seq:
       newseq += table[nucl]
    return "".join(newseq)


def compileMotifs(motifs, ignoreCase=True):
    """ input: list of motifs as dict motif->IUPAC, 
    returns: tuple -- list of regex objects, list of strings, list of descriptions """
    reList = []
    strList = []
    descList = []
    for m, desc in motifs.iteritems():
        name = m
        m = m.upper()
        revCompM = revComp(m)

        m = resolveIupac(m)
        strList.append(m)
        #print name, "+", m
        if ignoreCase:
            regex = re.compile(m, re.IGNORECASE)
        else:
            regex = re.compile(m)
        revCompM = resolveIupac(revCompM)

        if m!=revCompM:
            revCompRegex = re.compile(revCompM)
        else:
            revCompRegex = None

        reList.append((name,regex,revCompRegex))
        if desc==None:
            desc=m
        descList.append(desc)
    return reList, strList, descList


def scan(seq, regExList):
    """ input: sequence as string + list of tuples (motif name, regex, rev regex);
    output: dict iupac motif -> list of (position,strand) """
    assert(seq.find("-")==-1)
    assert(seq!=None)

    matches = {}
    for (motif, regEx, revCompRegEx) in regExList:
        posList = []
        findRegex(regEx, seq, posList, "+")
        if revCompRegEx != None:
            findRegex(revCompRegEx, seq, posList, "-")
        if len(posList)!=0:
            matches[motif]=posList

    return matches

def findMotifs(seq, motif):
    """ given seq as string + list of motifs as IUPAC strincs, 
        return occurences in format (start, end, strand, motif)
    """
    motifDict = {motif:motif}
    regExList, motRegexStr, descList = compileMotifs(motifDict)
    posList = []
    for (motif, regEx, revCompRegEx) in regExList:
        findRegex(regEx, seq, posList, "+")
        findRegex(revCompRegEx, seq, posList, "-")
    coordList = [(x,x+len(motif),strand, motif) for (x,strand) in posList]
    return coordList
