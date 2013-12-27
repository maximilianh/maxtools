#!/usr/bin/env python
from optparse import OptionParser
# python 2.4 default libraries
import sys, re, glob, urllib2, socket, os, doctest, logging, subprocess, os.path, fcntl, time
logger = logging.getLogger("t2gExtract")

# text2genome libraries
import t2gConfig, t2gConvert, namedtuple, maxTables, maxbio

# This program is: Copyright by Maximilian Haeussler maximilianh@gmail.com

# ======= CONSTANTS =======

# matches for these are removed from the file (=replaced by spaces)
xmlTagsRe  = re.compile('<.*?>')     # an xml tag
mathTypeRe = re.compile('MathType@[^ ]*') # a mathtype formula

# for cleaning/splitting the text files into words
nonLetterRe= re.compile(r'[\W]') # any non-alphanumeric character
digitRe    = re.compile('[0-9]')  # any digit
wordRe     = re.compile('[a-zA-Z]+') # any word

# for words consisting only of nucleodies
MINNUCLLEN = 3 # minimum length to add to stack
MINTRAILNUCLLEN = 17 # if string on stack is longer than this, MINNUCLLEN is ignored

# for non-nucleotide words
nonNucl       = re.compile("[^ACTGUactgu]") # regular expression that describes non-nucleotide letters (used for cleaning)
nuclRegex     = re.compile("[ACTGUactgu]")# regular expression that describes nucleotide letters
MINWORDLEN    = int(t2gConfig.get("text","MinWordLen", 19)) # minimum length of word to be processed
MINDNACONTENT = float(t2gConfig.get("text","MinDnaContent", 0.4)) # minimum content of nucleotide letters within a word

# for final filtering 
# minimum length of dna string to be output
MINDNALEN  = MINWORDLEN  
# maximum number of sequences per file to be output
MAXSEQS    = int(t2gConfig.get("text","maxSeqsPerDocument", 10000000)) 
# the minimum number of different letters in a word, to skip stuff like AGAGAG
MINDIFFLETTERS = 3 

CONTEXTLEN = 50 # characters around match that are output into tab sep file

# --- main wrapper called from command line tools

class MysqlGetter:
    """ get fulltext article data from mysql table """
    def __init__(self, connString, query):
        logger.info("Reading documents from mysql database %s with query %s" % (connString, query))
        self.table   = maxTables.SqlTableReader_BigTable(query, connString=connString)
        self.idCache = {} # status information for xml parser, need to know which shortId's have been used

    def generateDocuments(self):
        docCount = 0
        for row in self.table.generateRows():
            docCount+=1
            if docCount % 10000 == 0:
                logger.info("%d documents processed" % docCount)
            logger.debug("DocID %s" % row.documentId)

            if row.text!=None:
                logger.debug("using OCR text for extraction")
                text = row.text
            else:
                logger.debug("using XML text for extraction")
                text = row.xml

            yield row.documentId, text, row.xml

class FileGetter:
    """ get fulltext article data from text files in path """
    def __init__(self, folder):

        def findSubdirFiles(baseDir, extensions):
            """ Generator: traverse a baseDir and all subdirectories to find all files with certain extensions, extension is dot plus the extension, like ".xml"  """
            for root, dirs, files in os.walk(baseDir):
                for f in files:
                    if os.path.splitext(f)[1] in extensions:
                        yield os.path.join(root, f)

        if not os.path.isdir(folder):
            raise IOError("%s does not exist" % path)

        elif os.path.isfile(folder):
            logger.info("Only one file specified: Testing-mode, running on only one basename")
            filenames = [path]

        else:
            # get list of all basenames in this dir
            filenames = list(findSubdirFiles(folder, ['.xml','.txt']))

        baseNames = set([os.path.splitext(fname)[0] for fname in filenames])

        logger.info("Found %d basenames and %d files in directory %s\n" % (len(baseNames), len(filenames), folder))

        if len(baseNames)==0:
            logger.info("could not find any .xml or .txt files in directory %s\n" % folder)
            return

        self.baseNames = baseNames

    def generateDocuments(self):
        """ generator, yields tuple (docId, documentText, xml) """
        for fnameBase in self.baseNames:
            docId = filenameToId(fnameBase)

            # try several extensions, in this order
            txtnames = [fnameBase+".txt", fnameBase+".pdf.txt", fnameBase+".xml"]
            # find the first file with text
            text = None
            for txtname in txtnames:
                logger.debug("Trying %s" % txtname)
                occurrences = []
                if os.path.exists(txtname) and os.path.getsize(txtname)>=10:
                    logger.debug("Extracting from %s" % txtname)
                    text = open(txtname).read()
                    break
            
            if text==None: # skip doc if no text found
                logger.warn("documentId %s, could not find any text for this document" % fnameBase)
                continue

            xml = fnameBase+".xml"
            yield docId, text, xml

def extractPeptides(docGetter, tsvFile, fastaFile, peptideDetectCmdLine):
    """ iterates over files/directories yielded by a DocumentGetter and writes resulting peptides into tsvFile, fastaFile """

    docCount=0

    # start the perl script which detects protein sequences and return as a subprocess object

    for docId, text, xml in docGetter.generateDocuments():
        if docCount % 100 == 0:
            logger.info("%d documents processed" % docCount)
        docCount+=1

        text = cleanText(text)

        text.replace("-","")
        text.replace("~","")
        text.replace("(","")
        text.replace(")","")

        logger.debug("Executing %s" % peptideDetectCmdLine)
        proc = subprocess.Popen(peptideDetectCmdLine, bufsize=1000000, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        logger.debug("Piping article into peptide detector")
        stdoutData, stdinData = proc.communicate(text)

        pepLines = stdoutData.split("\n")
        logger.debug("Got: "+str(pepLines))

        for i in range(0, len(pepLines)):
            pepLine = pepLines[i]
            if len(pepLine)<=1:
                continue
            score, origText, pepSeq = pepLine.split("\t")
            if pepSeq.lower()=="elsevier":
                continue
            logger.debug("Found sequence: %s, original text %s, score %s" % (pepSeq, origText, score))

            data = [str(docId), "mysql", str(i), pepSeq, origText, score]
            tsvFile.write("\t".join(data)+"\n")
            logger.debug("Writing: "+str(data))
            fastaFile.write(">%s|%s\n%s\n" % (docId, i, pepSeq))
            #tsvFile.flush()
            #fastaFile.flush()

def extractNucleotides_files(path, tsvFile, fastaFile, metaInfoFile):
    """ iterate over files/directories (depending on path) and write results to tsv file and optionally fastaFile """
    #logger.debug("Start processing %s, writing output to %s and %s" % (path, tsvFilename, fastaFilename))

    idCache = {} # status information for xml parser
    if os.path.isdir(path):
        for docSeq, text in generateDocumentSequences(path):
            writeDocSeqToTsv ( tsvFile, docSeq)
            writeToFasta     ( fastaFile, docSeq)

            if metaInfoFile:
                assert(False) # PLEASE NOTE: PMC metainformation is broken for file extraction! have to return xml/text separately. remove metainfo in config and restart
                metaInfo = t2gConvert.parsePMCXmlAbstractInfo(docSeq.docId, text, idCache)
                maxbio.writeToTsv( metaInfoFile, metaInfo)

    elif os.path.isfile(path):
        logger.info("Only one file specified: Testing-mode, only outputting to stdout")
        occs   = nucleotideOccurrences(open(path).read())
        docSeq = DocumentSequences(filenameToId(path), path, occs)
        writeDocSeqToTsv (sys.stdout, docSeq)
    else:
        raise IOError("%s does not exist" % path)

def extractNucleotides_sql(connString, query, tsvFile, fastaFile, metaInfoFile):
    logger.info("Reading documents from mysql database %s with query %s" % (connString, query))
    table   = maxTables.SqlTableReader_BigTable(query, connString=connString)
    idCache = {} # status information for xml parser, need to know which shortId's have been used

    docCount = 0
    for row in table.generateRows():
        docCount+=1
        if docCount % 10000 == 0:
            logger.info("%d documents processed" % docCount)
        logger.debug("DocID %s" % row.documentId)

        if row.text:
            logger.debug("using OCR text for extraction")
            occs   = nucleotideOccurrences(row.text)
        else:
            logger.debug("using XML text for extraction")
            occs   = nucleotideOccurrences(row.xml)

        docSeq = DocumentSequences(int(row.documentId), "mysql", occs)
        writeToFasta(fastaFile, docSeq)
        writeDocSeqToTsv(tsvFile, docSeq)

        # convert xml to metaInfo
        if metaInfoFile:
            metaInfo = t2gConvert.parsePMCXmlAbstractInfo(row.documentId, row.xml, idCache)
            maxbio.writeToTsv( metaInfoFile, metaInfo)

# --- resolve filenames/locations -> document id --- 
class doiResolver:
# XX not used currently, but will be important for nature/elsevier
    def __init__(self):
        self.cache={}

    def convert(self, doi):
        if doi in self.cache:
            return self.cache[doi]

        url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?term=%s&email=maximilianh@gmail.com" % doi
        try:
            xml = urllib2.urlopen(url)
        except:
            return "PmidNotFound/HTTPError"

        for l in xml:
            if l.find("<Id>")!=-1:
                # <Id>16027735</Id>
                pmid = l.strip().replace("<Id>","").replace("</Id>", "")
                # strip off part after first _!
                self.cache[doi]=pmid
                return pmid

def filenameToId(fname):
    """ resolves a filename to a document id, returns an int """
    baseName = os.path.basename(fname)
    pmcId = int(os.path.splitext(baseName)[0].replace("PMC",""))
    return pmcId

# ==== SEQUENCE EXTRACTION =====

#  ... from directory

DocumentSequences = namedtuple.namedtuple("DocumentSequences", "docId, location, nucleotideOccurrences")

def generateDocumentSequences(folder):
    """ An iterator. Parses a folder with .txt and .xml files, but prefers txt files. Yields DocumentSequences"""

    # little helper
    def findSubdirFiles(baseDir, extensions):
        """ Generator: traverse a baseDir and all subdirectories to find all files with certain extensions, extension is dot plus the extension, like ".xml"  """
        for root, dirs, files in os.walk(baseDir):
            for f in files:
                if os.path.splitext(f)[1] in extensions:
                    yield os.path.join(root, f)

    # get list of all basenames in this dir
    filenames = list(findSubdirFiles(folder, ['.xml','.txt']))
    baseNames = set([os.path.splitext(fname)[0] for fname in filenames])

    logger.info("Found %d basenames and %d files in directory %s\n" % (len(baseNames), len(filenames), folder))

    if len(baseNames)==0:
        logger.info("could not find any .xml or .txt files in directory %s\n" % folder)
        return

    for fnameBase in baseNames:
        pmcId = filenameToId(fnameBase)

        # try several extensions, in this order
        txtnames = [fnameBase+".txt", fnameBase+".pdf.txt", fnameBase+".xml"]
        for txtname in txtnames:
            logger.debug("Trying %s" % txtname)
            occurrences = []
            if os.path.exists(txtname) and os.path.getsize(txtname)!=0:
                logger.debug("Extracting from %s" % txtname)
                text = open(txtname).read()
                occurrences = nucleotideOccurrences(text)
                logger.debug("Got %s" % str(occurrences))
                if len(occurrences)!=0:
                    break
        if len(occurrences)!=0:
            yield DocumentSequences(pmcId, txtname, occurrences), text

def writeDocSeqToTsv(tsvFile, docSeq):
    if tsvFile:
        logger.debug("Writing lines to tsvFile %s" % str(tsvFile))
        for occ in docSeq.nucleotideOccurrences:
            data = [docSeq.docId, docSeq.location, occ.seqId, occ.seq, occ.start, occ.end, occ.context, occ.partCount, int(occ.tainted)]
            data = [str(x) for x in data]
            tsvFile.write("\t".join(data)+"\n")

def writeToFasta(faFile, docSeq):
    logger.debug("Writing lines to fastaFile %s" % str(faFile))
    for occ in docSeq.nucleotideOccurrences:
        faFile.write(">%d|%d\n%s\n" % (docSeq.docId, occ.seqId, occ.seq))

# ... from file

NucleotideOccurrence = namedtuple.namedtuple("NucleotideOccurrence", "start, end, seq, seqId, context, partCount, tainted")

def replaceWithSpaces(regex, string):
    """ replaces all occurrences of regex in string with spaces 
    >>> replaceWithSpaces(xmlTagsRe, "<test> nonTag <test>")
    '       nonTag       '
    """
    def toSpaces(matchObject):
        return "".join([" "]*(matchObject.end(0) - matchObject.start(0)))
    return regex.sub(toSpaces, string) 

class MatchStack:
    """ a stack of matches to nucleotide regex objects in a text """
    def __init__(self):
        self.reset()

    def cleanPush(self, match):
        """ clean a word from non-nucleotides, then push onto stack"""
        self.tainted=True
        self.push(match, clean=True)

    def push(self, match, clean=False):
        newString = match.group(0)

        if clean:
            newString = nonNucl.sub("", newString) # remove all non-nucleotide characters

        newString = newString.replace("u", "t").replace("U", "T") # RNA -> DNA
        logger.log(5, "stack push: %s" % newString)
        self.strings.append(newString)
        self.charLen += len(newString)

        self.start = min(self.start, match.start())
        self.end   = max(self.end, match.end())

    def diffLetters(self):
        longString = "".join(self.strings)
        longString = longString.lower()
        return len(set(longString))
        
    def hasLongSequence(self):
        return (self.charLen >= MINDNALEN and self.diffLetters() >= MINDIFFLETTERS)

    def getOcc(self, seqId, text):
        nuclString         = "".join(self.strings)
        logger.log(5, "generating sequence from stack: id %d, sequence %s" % (seqId, nuclString))

        validStart = max(0, self.start-CONTEXTLEN)
        validEnd   = min(len(text), self.end+CONTEXTLEN)
        context    = text[validStart:validEnd]
        partCount  = len(self.strings)
        return NucleotideOccurrence(self.start, self.end, nuclString, seqId, context, partCount, self.tainted)

    def reset(self):
        logger.log(5, "stack reset")
        self.strings = []
        self.start   = 999999999 # start pos of first element
        self.end     = 0         # end position of last element
        self.charLen = 0         # length of whole stack in characters
        self.tainted = False     # has ASCII-cleaning been used to populate any element in the stack?


def cleanText(text):
    # clean: xml tags and mathtype -> spaces
    cleanText = replaceWithSpaces(xmlTagsRe, text)
    cleanText = replaceWithSpaces(mathTypeRe, text)
    # clean: non-letters -> spaces
    cleanText = nonLetterRe.sub(" ", cleanText)
    cleanText = digitRe.sub(" ", cleanText)
    return cleanText

def nucleotideOccurrences(text):
    """ Parse out all nucleotide-like strings from xml/ascii text and return them as a list of nucleotide occurence-records
    Example:

    >>> nucleotideOccurrences("test test caccatgacacactgacacatgtgtactgtg")[0]
    NucleotideOccurrence(seqId=0, seq='caccatgacacactgacacatgtgtactgtg', start=10, end=41, context='test test caccatgacacactgacacatgtgtactgtg', partCount=1, tainted=False)
    >>> nucleotideOccurrences("test test tga tga cac atg tgt act gtg a")[0].seq
    'tgatgacacatgtgtactgtga'
    >>> nucleotideOccurrences("bla bla actg ttt tcactybaactbacbatactbatcgactgactgactgtactcctacgatgcgtactacttacghhh")[0].seq
    'actgttttcactaactacatactatcgactgactgactgtactcctacgatgcgtactacttacg'

    """

    textClean = cleanText(text)

    stack = MatchStack()
    seqId = 0
    occurences = []

    for wordMatch in wordRe.finditer(textClean):
        word = wordMatch.group(0)
        logger.log(5, "word: %s" % word)

        dnaContent = float(len(nuclRegex.findall(word))) / len(word)

        if dnaContent==1.0 and (len(word)>=MINNUCLLEN or stack.charLen>MINTRAILNUCLLEN):
            # word is only nucleotides: either long enough or already enough
            # chars in stack
            stack.push(wordMatch)
            continue
        elif len(word) >= MINWORDLEN and dnaContent >= MINDNACONTENT:
            # long words with enough DNA in them
            stack.cleanPush(wordMatch)
        else:
            # anything else triggers a stack output and reset
            if stack.hasLongSequence():
                occurences.append(stack.getOcc(seqId, textClean))
                seqId+=1
            stack.reset()

    # if document finishes with a non-empty stack: empty it
    if stack.hasLongSequence():
        occurences.append(stack.getOcc(seqId, textClean))

    if len(occurences) > MAXSEQS:
        logger.log(5, "too many sequences in paper, skipping whole document")
        return []
    else:
        return occurences

def pubtools_extractNucleotides(metaInfo, text):
    """ interface for pubtools """
    result = []
    for line in nucleotideOccurrences(text):
        line = [str(x) for x in line]
        line = "\t".join(line)
        result.append(line)
    return "\n".join(result)

# ----- 
if __name__ == "__main__":
    import doctest
    doctest.testmod()

