from sys import *
from optparse import OptionParser
import logging, os, cgi
import tabfile, pmcArticle, t2gConfig, util, maxTables, maxbio, unicodeConvert

# for pmc xml parser
from xml.sax import make_parser
from xml.sax.handler import ContentHandler
import xml.sax.handler, os.path
import glob, re
import util
from StringIO import StringIO

logger = logging.getLogger("t2gConvert")

try:
    import unidecode
except:
    logging.debug("No unidecode installed")

try:
    import MySQLdb
except:
    #logger.warn("MySQLdb module is not installed on this computer, you cannot read articles from mysql databases")
    pass

def featuresToBed(featureFile, taxonToUcsc, pmcToKeyFile, outputDir, sqlConnDict, forcePosList, asCustomTrack, writeDefaultPos, convertChroms):
    """ convert features.tsv from t2gProcess.process to bed format for the UCSC genome browser """

    # this is ugly but so are genomes: Maps from Ensembl to UCSC chromosome names
    chromMap = {
            "chrMtDNA" : "chrM", # ce2
            "chrUextra" : "chrU", # ce2
            "E22C19W28_E50C23" : "chrE22C19W28_E50C23", #galgal2
            "chrMT" : "chrM",  # mm9
            "X" : "chrX",
            "Y" : "chrY", 
            "Uextra" : "chrUextra",
            "U" : "chrU",
            "dmel_mitochondrion_genome" : "chrM", # dm2
            }

    latinToArabic = {
            "chrI" : "chr1", 
            "chrII" : "chr2", 
            "chrIII" : "chr3", 
            "chrIV" : "chr4", 
            "chrV" : "chr5", 
            "chrVI" : "chr6", 
            "chrVII" : "chr7", 
            "chrVIII" : "chr8", 
            "chrIX" : "chr9", 
            "chrX" : "chr10", 
            "chrXI" : "chr11", 
            "chrXII" : "chr12", 
            "chrXIII" : "chr13", 
            "chrXIV" : "chr14" ,
            "chrXV" : "chr15" ,
            "chrXVI" : "chr16" ,
            "chrMito" : "chrM" 
            }

    arabicToLatin = {
            "I"    : "chr1",
            "II"   : "chr2",
            "III"  : "chr3",
            "IV"   : "chr4",
            "V"    : "chr5",
            "VI"   : "chr6",
            "VII"  : "chr7",
            "VIII" : "chr8",
            "IX"   : "chr9",
            "X"    : "chr10",
            "XI"   : "chr11",
            "XII"  : "chr12",
            "XIII" : "chr13",
            "XIV"  : "chr14" ,
            "XV"   : "chr15" ,
            "XVI"  : "chr16" ,
            "MT"   : "chrM",
            "Mito" : "chrM"
            }

    def ensemblToUcsc(chrom):
        # convert ensembl to UCSC chromosomes
        # a very ugly procedure but not our fault
        # (Ensembl itself doesn't care, they don't convert stuff like HSCHR4, creating invalid links)
        # (UCSC doesn't care either, links to Ensembl often don't work when not on standard chromosomes)
        # I spent two hours fixing most of these links, a very ugly work
        
        if ucscDb.startswith("sacCer1") and chrom=="2-micron":
            return None
        elif ucscDb.startswith("sacCer1"):
            chrom = chrom.replace("chr","")
            chrom = arabicToLatin[chrom]
        elif ucscDb.startswith("dm") and chrom.startswith("U"):
            return None
        elif ucscDb.startswith("mm") and chrom.startswith("NT_"):
            return None
        elif ucscDb.startswith("hg") and "_" in chrom:
            return None
        elif ucscDb.startswith("galGal") and "_" in chrom:
            return None
        elif ucscDb.startswith("hg") and chrom.startswith("GL"):
            return None
        elif ucscDb.startswith("rheMac") and len(chrom)>10:
            return None
        elif ucscDb.startswith("tetNig") and ("_" in chrom or "E64" in chrom):
            return None
        elif ucscDb.startswith("galGal") and "_" in chrom:
            return None
        # on bosTau4
        elif ucscDb.startswith("bosTau") and chrom.startswith("AAFC"):
            return None

        # if the chrom is a simple number we have to add "chr"
        if chrom[0].isdigit():
            # some hacks to accommodate ensembl/UCSC pecularities
            if ucscDb=="ci2":
                chrom = "chr%02d" % int(chrom.replace("q","").replace("p","")) + chrom[-1]
            else:
                chrom = "chr"+chrom
        elif chrom[0] in ["I","V","M","X","Y","U","Z","W"]:
            chrom = "chr"+chrom


        if chrom in chromMap:
            chrom = chromMap[chrom]
        return chrom

        # END OF ENSEMBLTOUCSC

    # if pmcTOKeyfile looks like sql query, parse it into a dict
    if pmcToKeyFile==None:
        pmcToKey = {}
    elif pmcToKeyFile.lower().startswith("select"):
        pmcToKey = maxTables.SqlTableReader(pmcToKeyFile, **sqlConnDict).asDict()
    else:
        pmcToKey = tabfile.slurpdict(pmcToKeyFile, keyAsInt=True, comments=True)

    # reading input files and opening target files -> bedFiles
    #logger.info("Reading %s" % taxonToUcscFile)
    bedFiles = {}
    for taxon, ucscDb in taxonToUcsc.iteritems():
        # UCSC's fr2 has nothing in common with ensembl's fugu sequences
        if ucscDb=="fr2":
            continue
        bedname = os.path.join(outputDir, ucscDb+".bed")
        fh = open(bedname, "w")
        if asCustomTrack:
            fh.write('track name=Text2Genome description="Text2Genome Fulltext DNA" visibility=3 htmlUrl=http://www.text2genome.org url=%s?displayid=$$\n' % t2gConfig.detailsUrl)
        bedFiles[ucscDb] = fh

    logger.info("Reading features.tsv into memory")
    indexedFeatures = {}
    i = 0
    dropped =0
    for line in open(featureFile):
        if line.startswith("#"):
            continue
        line = line.strip()
        ft = line.split("\t")
        genomeId, chrom, start, end, pmcId, seqId, score, groupId = ft[:8]
        pmcId=int(pmcId)
        groupId=int(groupId)

        ft = ( genomeId, chrom, int(start), int(end), pmcId, int(seqId), int(float(score)), groupId )
        if genomeId in taxonToUcsc:
            indexedFeatures.setdefault(pmcId, {}).setdefault(groupId, []).append(ft)
            i+=1
        else: 
            dropped+=1
    logger.info("Read %d features, dropped %d features as genome is not mapped to UCSC db ig" % (i, dropped))

    # format of feature is: genomeId, chrom, start, end, pmcId, seqId, score, groupId = ft
    # this could be better realized with dictionaries but I think arrays are faster and the code
    # is still quite readable
    GENOME=0
    CHROM=1
    START=2
    END=3
    PMCID=4
    SEQID=5

    # sample feature for default positions
    sampleDefPos = {}

    logger.info("converting features and writing to files in %s" % outputDir)
    for pmcId, groupFeatures in indexedFeatures.iteritems():
        for groupId, features in groupFeatures.iteritems():

            # data check
            chroms = [f[CHROM] for f in features]
            chroms = set(chroms)
            assert(len(chroms)==1) # all features of a group have to be located on the same chromosome
            chrom = chroms.pop()

            genomeId = features[0][GENOME]
            ucscDb = taxonToUcsc[genomeId]

            if ucscDb not in bedFiles:
                dropped+=len(features)
                continue

            if convertChroms:
                chrom = ensemblToUcsc(chrom)
                if chrom==None:
                    continue

            features.sort(key=lambda f: f[START]) # sort fts by start pos

            # generate bed block sizes for bed feature line

            # generate bitmask of occupied positions
            minStart = min([f[START] for f in features])
            maxEnd = max([f[END] for f in features])
            mask = [0]*(maxEnd-minStart)
            for f in features:
                for i in range(f[START]-minStart,f[END]-minStart):
                    mask[i]=1
                #print f, mask

            blockStarts = []
            blockSizes = []
            # search for consec stretches of 1s
            lastStart=None
            wasOne=False
            wasZero=True
            #print mask
            for i in range(0, len(mask)):
                #print i
                if mask[i]==1 and wasZero:
                    #print "start"
                    blockStarts.append(i)
                    wasZero=False
                    wasOne=True
                    lastStart=i
                if mask[i]==0 and wasOne:
                    #print "end"
                    blockSizes.append(i-lastStart)
                    wasZero=True
                    wasOne=False
                    lastStart=None
            if lastStart!=None:
                blockSizes.append(len(mask)-lastStart)

            blockStarts = [str(x) for x in blockStarts]
            blockSizes = [str(x) for x in blockSizes]
            seqIds = set([f[SEQID] for f in features])
            seqIds = [str(x) for x in seqIds]
            
            name = str(pmcToKey.get(pmcId, pmcId))
            score = 0
            bedFeature = [chrom, str(minStart), str(maxEnd), name ,str(score), "+", str(minStart), str(maxEnd), "0", str(len(blockStarts)), ",".join(blockSizes), ",".join(blockStarts)]
            if not asCustomTrack:
                bedFeature.append(",".join(seqIds))

            bedFh  = bedFiles[ucscDb]
            bedFh.write("\t".join(bedFeature))
            bedFh.write("\n")
            sampleDefPos[ucscDb]="%s:%d-%d" % (chrom, minStart, maxEnd)

    logger.info("Created bed files for the following UCSC assemblies: "+", ".join(bedFiles.keys()))

    if forcePosList and writeDefaultPos:
        for db, pos in forcePosList:
            logger.debug("Manually setting default position for %s to %s" % (db, pos))
            sampleDefPos[db]=pos

    if sqlConnDict!=None and sqlConnDict!="" and len(sqlConnDict)!=0 and writeDefaultPos:
        logger.info("Writing default positions to mysql database, table ucscDbdb")
        db = MySQLdb.connect(**sqlConnDict) 
        cursor = db.cursor()

        for dbname, defaultPos in sampleDefPos.iteritems():
            cursor.execute('update ucscDbdb set defaultPos="%s" where name="%s";' % (defaultPos, dbname))

        cursor.close ()
        db.commit()
        logger.debug("Wrote the following default positions: %s" % (str(sampleDefPos)))

def convertBlastFiles(blastDirs, genomeToTaxFile, tempFilename, outFilename, fileFormat):
    """ collect all *.blast files from blastDirs and subdirs, map to taxon ids using genomToTaxFile and write results to outfh fileobject 
    fileFormat can be blat or blast """

    outfh               = maxbio.openFile(tempFilename, "w")

    # read genome -> taxid map 
    if genomeToTaxFile:
        orgToNum = {}
        logger.info( "Reading genome -> number table from %s" % genomeToTaxFile)
        logger.info( "Expected input fields are: (GenomeName, other field,..., other field, genomeId)")

        for l in open(genomeToTaxFile):
            if l.startswith("#"):
                continue
            fs = l.strip().split("\t")
            genome = fs[0]
            num = fs[-1]
            num = int(num)
            genome=genome.lower().replace(" ", "_")
            orgToNum[genome]=num
    else:
        logger.info( "No genome map specified, results are not genome-based")
        orgToNum=None

    dropFileCount=0
    lineCount=0
    finishedTaxons = set()
    lastDir = None

    # read blast files
    if fileFormat=="blast":
        ext = ".blast"
    elif fileFormat=="blat":
        ext = ".psl"
    elif fileFormat=="bwa":
        ext = ".sam"
    else:
        assert(False) # wrong file format parameter

    for blastDir in blastDirs:
        logger.info("Searching for files with extension %s in directory %s" % (ext, blastDir)) 
        files = list(util.findSubdirFiles(blastDir, ext))
        logger.info("Found %d files" % len(files))

        # convert blast files
        files = list (files)
        for fname in files:
            # convert organism to taxid, skip if not possible
            dirname = os.path.dirname(fname)
            org = os.path.basename(dirname)
            org = org.lower()

            if orgToNum!=None:
                if org in orgToNum:
                    orgNum = orgToNum[org]
                else:
                    # try to find any organism from genome list in filename
                    found = False
                    for dbOrg in orgToNum:
                        if dbOrg.replace(" ", "_").lower() in fname.lower():
                            orgNum = orgToNum[dbOrg]
                            logger.info("Found orgName %s in filename %s, using organism id %s" % (dbOrg, fname, str(orgNum)))
                            found=True
                            break

                    if not found:
                        logger.warn("warning: could not resolve filename %s to taxid, dropping this file (recognized organism %s)" % (fname, org))
                        dropFileCount+=1
                        continue
            else:
                orgNum=-1

            # check if not already processed AND in different directory (blast creates several indices per directory), skip if yes
            #if orgNum in finishedTaxons and dirname!=lastDir:
                #print("warning: already processed this taxon id %d, skipping input file %s)" % (orgNum, fname))
                #continue
            finishedTaxons.add(orgNum)
            lastDir = dirname

            # convert lines
            f = open(fname, "r")
            #print "Reading %s, writing hits to %s"%(fname,outFile)
            if orgNum!=-1:
                logger.info( "Reading %s, genomeID %d"%(fname, orgNum))
            else:
                logger.info( "Reading %s, not linked to any genome id" % (fname))

            if fileFormat=="bwa":
                tp = maxTables.TableParser(fileType="sam")

            for l in f:
                lineCount+=1
                # parse blast line
                fs = l.strip().split("\t")
                if fileFormat=="blast":
                    # example
                    # 11495631        chr1    100.00  23      0       0       1       23      25500772        25500750        2e-05   46.1
                    srcId, trgId, perc, length, dummy, dummy, dummy, length, trgStart, trgEnd, eVal, score = fs
                elif fileFormat=="blat":
                    # psl-format from http://genome.ucsc.edu/FAQ/FAQformat.html#format2
                    # matches - Number of bases that match that aren't repeats
                    # misMatches - Number of bases that don't match
                    # repMatches - Number of bases that match but are part of repeats
                    # nCount - Number of 'N' bases
                    # qNumInsert - Number of inserts in query
                    # qBaseInsert - Number of bases inserted in query
                    # tNumInsert - Number of inserts in target
                    # tBaseInsert - Number of bases inserted in target
                    # strand - '+' or '-' for query strand. For translated alignments, second '+'or '-' is for genomic strand
                    # qName - Query sequence name
                    # qSize - Query sequence size
                    # qStart - Alignment start position in query
                    # qEnd - Alignment end position in query
                    # tName - Target sequence name
                    # tSize - Target sequence size
                    # tStart - Alignment start position in target
                    # tEnd - Alignment end position in target
                    # blockCount - Number of blocks in the alignment (a block contains no gaps)
                    # blockSizes - Comma-separated list of sizes of each block
                    # qStarts - Comma-separated list of starting positions of each block in query
                    # tStarts - Comma-separated list of starting positions of each block in target
                    # 23      0       0       0       0       0       0       0       +       11075971|1      2299    2248    2271    scaffold_281    111378  17336   17359   1       23,     2248,   17336,
                    matches, misMatches, repMatches, nCount, qNumInsert, qBaseInsert, tNumInsert, tBaseInsert, strand, qName, qSize, qStart, qEnd, tName, tSize, tStart, tEnd, blockCount, blockSizes, qStarts, tStarts = fs

                    score = matches
                    perc = "%2.2f" % ((float(matches)+float(repMatches)) / (float(matches)+float(misMatches)+float(repMatches)) * 100.0)

                    trgId = tName
                    srcId = qName
                    trgStart = tStart
                    trgEnd   = tEnd

                elif fileFormat=="bwa":
                    if fs[0].startswith("@"):
                        continue
                    if len(fs)==11:
                        fs.append("")
                    row = tp.parseTuple(fs)
                    tuple = maxTables.samToBed(row)
                    if tuple==None:
                        continue

                    chrom, start, end, name, score, strand = tuple
                    srcId = name
                    trgId = chrom
                    perc = "0"
                    length = start-end
                    trgStart, trgEnd=start, end
                    
                else:
                    assert(False) # file format not found?


                trgEnd = int(trgEnd)
                trgStart = int(trgStart)
                if trgEnd < trgStart:
                    trgStart, trgEnd = trgEnd, trgStart


                fs = srcId.split("|")
                srcId = fs[0]
                srcSeq = fs[1]
                pmcId = srcId.replace("PMC", "")
                pmcId = pmcId.replace(".txt", "")
                pmcId = pmcId.replace(".pdf", "")
                data = [pmcId, str(orgNum), srcSeq, trgId, str(trgStart), str(trgEnd), str(score), str(perc)]
                outfh.write("\t".join(data)+"\n")

    outfh.close()

    logger.info("BlastHit output table format is (pmcId, genomeId, seqId, chrom, start, end, score, percentId)")
    logger.info("blastConvert : blast files dropped because of unresolvable species name %d, filesDropped=%d" % (dropFileCount, dropFileCount))
    logger.info("blastConvert : processed %d blast matches, blastMatches=%d" % (lineCount, lineCount))
    logger.info("Now sorting the file with the UNIX sort command")

    cmdLine = "sort -n %s | uniq > %s" % (tempFilename, outFilename)
    logger.info(cmdLine)
    ret = os.system(cmdLine)

    if ret==0:
        logger.info("Sorting finished, no error")
    else:
        logger.info("Error occured while sorting")

# ----- xml parsing of articles
nonLetterRe = re.compile('\W', re.U) # non-alphanumeric character in unicode set, uff.

class ArticleHandler(ContentHandler):
    def __init__ (self, uniqueIdCounts): 
        self.uniqueIdCounts = uniqueIdCounts
        self.surnames = []
        self.givenNames = []
        self.isJournalTitle = False
        self.isJournalId = False
        self.isPmid = False
        self.isSubject = False
        self.isArticleTitle= False
        self.isSurname = False
        self.isGivenNames = False
        self.isPubDate=False
        self.isMonth=False
        self.isYear=False
        self.isVolume=False
        self.isIssue=False
        self.isFpage=False
        self.isLpage=False
        self.isFront=False
        self.isAbstract=False
        self.isArticle=False

        self.subject = "";
        self.journalTitle = ""
        self.journalId = ""
        self.pmid=""
        self.abstract=""
        self.articleTitle = "";
        self.month = ""
        self.lpage = ""
        self.years = []
        self.volume = ""
        self.issue = ""
        self.fpage = ""

    def startElement(self, name, attrs):
       if name == 'article':
           self.isArticle=True

       if name == 'front':
           self.isFront = True;

       if self.isFront and name=='abstract':
            self.isAbstract = True

       if self.isAbstract:
           if name == 'title':
               self.abstract += '<b>'

       if self.isFront:
           if name == 'journal-title':     
             self.isJournalTitle = True
           elif name == 'journal-id' and attrs.get("journal-id-type", "")=="nlm-ta":     
             self.isJournalId = True
           elif name == 'article-id':
             if attrs.get("pub-id-type", "")=="pmid":
                self.isPmid=True
           elif name == 'subject':
             self.isSubject = True;
           elif name == 'article-title':
             self.isArticleTitle = True;
           elif name == 'contrib' and attrs.get("contrib-type", "")=="author":
             self.isContrib = True;
           elif name == 'surname':
             self.isSurname = True;
             self.surnames.append("")
           elif name == 'given-names':
             self.isGivenNames = True;
             self.givenNames.append("")
           elif name == 'volume':
             self.isVolume = True;
           elif name == 'issue':
             self.isIssue = True;
           elif name == 'fpage':
             self.isFpage = True;
           elif name == 'lpage':
             self.isLpage = True;
           elif name == 'pub-date' and (attrs.get("pub-type","")=="ppub" or attrs.get("pub-type", "")=="epub"):
             self.isPubDate = True
           elif self.isPubDate:
               if name == 'month':
                 self.isMonth = True;
               elif name == 'year':
                 self.isYear = True;

    def endElement(self, name):
       if name == 'journal-id':
         self.isJournalId = False
       if name == 'abstract' and self.isFront:
         self.isAbstract = False
       if self.isAbstract and name=="title":
         self.abstract += "</b> "
       if name == 'journal-title':
         self.isJournalTitle = False
       if name == 'article-id':
         self.isPmid = False
       if name == 'subjet':
         self.isSubject = False
       if name == 'front':
         self.isFront = False
       if name == 'article-title':
         self.isArticleTitle = False
       if name == 'contrib':
         self.isContrib = False
       if name == 'surname':
         self.isSurname = False
       if name == 'given-names':
         self.isGivenNames = False
       if name == 'pub-date':
         self.isPubDate = False
       if name == 'month':
         self.isMonth = False
       if name== 'year':
         self.isYear = False
       if name=='volume':
           self.isVolume = False
       if name=='issue':
           self.isIssue = False
       if name=='fpage':
           self.isFpage = False
       if name=='lpage':
           self.isLpage = False

    def characters (self, ch):
        if self.isFront:
           ch = ch.strip("\n")
           ch = ch.replace("\t"," ")
           if     self.isJournalTitle:
                  self.journalTitle += ch
           if     self.isAbstract:
                  self.abstract += ch
           if     self.isJournalId:
                  self.journalId += ch
           if     self.isPmid:
                  self.pmid = str(int(ch)) # strip leading 000 from pmid
           if     self.isSubject:
                  self.subject += ch
           if     self.isArticleTitle:
                  self.articleTitle += ch
           if     self.isSurname:
                  self.surnames[-1]+=ch
           if     self.isGivenNames:
                  self.givenNames[-1]+=ch
           if     self.isPubDate:
               if self.isMonth:
                  self.month = ch
               if self.isYear:
                  self.years.append(ch)
           if     self.isVolume:
                  self.volume += ch
           if     self.isIssue:
                  self.issue += ch
           if     self.isFpage:
                  self.fpage += ch
           if     self.isLpage:
                  self.lpage += ch

    def toString(self, pmcId):
        authors = []
        for i in range(0, min(len(self.surnames), len(self.givenNames))):
            authors.append(self.surnames[i]+", "+self.givenNames[i])
        authorsStr = "|".join(authors)

        if "journalId" in self.__dict__:
            journal = self.journalId
        else:
            journal = self.journalTitle

        if "volume" not in self.__dict__:    
            pages = ""
        elif "issue" not in self.__dict__:
            pages = "%s" % self.volume
        elif "fpage" not in self.__dict__:
            pages = "%s(%s)" % (self.volume, self.issue)
        elif "lpage" not in self.__dict__:
            pages = "%s(%s):%s" % (self.volume, self.issue, self.fpage)
        else:
            pages = "%s(%s):%s-%s" % (self.volume, self.issue, self.fpage, self.lpage)

        years = [int(y) for y in self.years if y!='']  # try to find year when article was first published
        #stderr.write(pmcId+years)
        if len(years)>0:
            minYear = str(min(years))
        else:
            minYear=0

        # create a unique id which does not exist yet
        if len(authors)==0:
            uniqueId = pmcId
        else:
            if len(self.surnames)>0:
                firstAuthor = self.surnames[0]
                firstAuthor = nonLetterRe.sub("", firstAuthor) # remove nonalphanumeric characters
                #firstAuthor = util.remove_accents(firstAuthor) # remove accents
                firstAuthor = unidecode.unidecode(firstAuthor) # remove accents
                # test if this works with pmcid 1936429
            else:
                firstAuthor = "NoAuthor"
            stem = firstAuthor+minYear
            stemCount = self.uniqueIdCounts.get(stem, 0)
            if stemCount == 0:
                self.uniqueIdCounts[stem]=1
                suffix=""
            else:
                count = self.uniqueIdCounts[stem]
                self.uniqueIdCounts[stem]+=1 
                suffix = util.baseN(count) 
            uniqueId = stem+suffix

        # print result
        line = [pmcId, uniqueId, authorsStr, self.articleTitle, journal, minYear, pages, self.pmid, self.abstract]
        return line

def parsePMCXmlAbstractInfo(pmcId, string, idCache):
    artHandler = ArticleHandler(idCache)

    parser = make_parser()   
    parser.setFeature(xml.sax.handler.feature_validation, 0) 
    parser.setFeature(xml.sax.handler.feature_namespaces, 0) 
    parser.setFeature(xml.sax.handler.feature_external_pes, 0) 
    parser.setFeature(xml.sax.handler.feature_external_ges, 0) 

    parser.setContentHandler(artHandler)
    try:
        parser.parse(StringIO(string))
        if artHandler.isArticle:
            return artHandler.toString(pmcId)
    except :
        logger.info("error: exception %s, details %s, PMCID %s\n" % (exc_info()[0],exc_info()[1], pmcId))

def exportToUcscTables(dbDict, taxonToUcsc, outputDir):
    logger.info("If this process hangs, try to restart mysqld")
    for taxon, db in taxonToUcsc.iteritems():
        logger.info("Creating UCSC tables for taxon %s, ucsc db %s" % (str(taxon), str(db)))
        #SET NAMES utf8; SET CHARACTER SET utf8; SET character_set_connection=utf8; 
        articleQuery = 'SELECT distinct pmcId as docId, displayId, replace(authors,"|","; ") as authors, title, journal, year, pages, pmid, abstract from refs JOIN features USING (pmcId) WHERE features.genomeId=\'%s\';' % taxon
        seqQuery = "select distinct pmcId as docId, sequences.seqId, sequence from sequences JOIN features USING (pmcId) JOIN refs USING (pmcId) WHERE features.genomeId='%s'" % taxon
        logger.debug("getting articles...")
        artTr = maxTables.SqlTableReader_BigTable(articleQuery, **dbDict)
        logger.debug("getting sequences...")
        seqTr = maxTables.SqlTableReader_BigTable(seqQuery, **dbDict)

        artName = os.path.join(outputDir, "t2gArticle",db+".tab")
        artFile = open(artName, "w")

        logger.debug("outputting article table...")
        for row in list(artTr.generateRows()):
            #text = util.removeAccents(text)
            text = "\t".join([str(f) for f in row]) +"\n"
            try:
                text = unicode(text, "utf8")
            except:
                text = unicode(text, "latin1")
            text = text.replace("<b>Images:</b>","")
            text = text.replace("<b>","<br><b>").replace("</b>",":</b> ")
            #text = text.encode('latin1', 'xmlcharrefreplace')
            text = unicodeConvert.string_to_ncr(text)
            artFile.write(text)

        seqName = os.path.join(outputDir, "t2gSequence",db+".tab")
        seqFile = open(seqName, "w")

        logger.debug("outputting sequences table...")
        for row in seqTr.generateRows():
            text = "\t".join([str(f) for f in row])+"\n"
            seqFile.write(text)

def tabToFasta(inFilename, outFilename, linnaeusFilename, genomesFilename, numArtIdField=0, seqIdField=1, artIdField=2, seqField=5, maxTotalSeqLen=1000000, maxSeqLen=100000, maxSeqCount=100):
    """ convert tab-sep file from fttools to fasta file, linnaeusFile has docId<tab>taxonId format and genomesFilename has <genomeId> format """

    if linnaeusFilename!=None:
        logging.debug("Reading LINNAEUS files")
        docToTaxon      = tabfile.slurpdictset(linnaeusFilename)
        logging.debug("Reading sequenced genome list")
        sequencedTaxons = set(tabfile.slurplist(genomesFilename))
        taxonFileObjCache = {}

        logging.debug("Filtering LINNAEUS data by sequenced genome list")
        # filter the taxonIds by sequenced genome Ids
        filteredDocToTaxon = {}
        allSequencedTaxonIds = set()

        for docId, taxonIds in docToTaxon.iteritems():
            filtTaxons = set()
            for taxonId in taxonIds:
                if taxonId in sequencedTaxons:
                    filtTaxons.add(taxonId)
                    allSequencedTaxonIds.add(taxonId)
            filteredDocToTaxon[docId]=filtTaxons

        docToTaxon = {} # no need to waste memory with this
        logging.debug("Got data for %d sequenced genomes" % len(allSequencedTaxonIds))

        # open filehandles for these
        taxonToFile = {}
        for taxonId in allSequencedTaxonIds:
            taxonToFile[taxonId] = open(outFilename+"."+str(taxonId)+".fa", "w")

    else:
        docToTaxon=None
        outFileObj = open(outFilename, "w")

    colCount = len(open(inFilename).readline().split("\t"))

    br = maxTables.BlockReader(open(inFilename), 0, mustBeSorted=False)

    for artId, block in br.readNext():
        sequences = []
        seqSet = set()
        internalId = None
        externalId = None
        for fs in block:
            internalId = fs[numArtIdField]
            externalId = fs[artIdField]
            seq = fs[seqField]
            seqId = fs[seqIdField]

            letterCount = len(set(seq))
            if letterCount<=2:
                logging.debug("Skipping sequence %s of article %s, not more than 2 letters in string: %s" % (seqId, externalId, seq))
            elif seq in seqSet:
                logging.debug("Skipping sequence %s of article %s, already seen before: %s" % (seqId, externalId, seq))
            elif len(seq)>maxSeqLen:
                logging.debug("Skipping sequence %s of article %s; longer than %d bp: %s" % (seqId, externalId, maxSeqLen, seq))
            else:
                sequences.append ( (seqId, seq) )
                seqSet.add(seq)

        seqCount = len(sequences)
        totalSeqLen = sum([len(y) for x,y in sequences])

        if totalSeqLen > maxTotalSeqLen:
            logging.debug("Skipping article %s, total sequence length %d > %d" % (externalId, totalSeqLen, maxTotalSeqLen))
        elif seqCount>maxSeqCount:
            logging.debug("Skipping article %s, total sequence count %d > %d" % (externalId, seqCount, maxSeqCount))
        else:
            for seqId, seq in sequences:
                if docToTaxon==None:
                    outFileObj.write(">%s|%s\n%s\n" % (internalId, seqId, seq))
                else:
                    taxonIds = filteredDocToTaxon.get(internalId, None)
                    #logging.debug("TaxonIds are %s" % str(taxonIds))
                    if taxonIds==None:
                        taxonIds = allSequencedTaxonIds
                    for taxonId in taxonIds:
                        outFileObj = taxonToFile[taxonId]
                        outFileObj.write(">%s|%s\n%s\n" % (internalId, seqId, seq))



snpRegex = re.compile('rs[0-9]{5,8}')
space = re.compile(r'\W+')

def pubtools_snpSearch(metaInfo, text):
    matches = snpRegex.finditer(text)
    lines = []
    textLen = len(text)
    for match in matches:
        start = match.start()
        end = match.end()
        matchStr = text[start:end]

        contextStart = max(0, start-100)
        contextEnd  = min(end+100, textLen)
        context = text[contextStart: contextEnd].replace("\n"," ")
        context = space.sub(" ", context) # collapse spaces
        fields = [str(start), str(end), matchStr, context]
        line = "\t".join(fields)
        lines.append(line)
    return "\n".join(lines)
    
rtPcrRegex = re.compile('(qpcr)|(q-pcr)|(qrt-pcr)|(quantitative pcr)|(quantitative poly.*)|(quantitative polymerase)|(quantitative real-time)|(realtime pcr)|(reverse transcription polymerase)|(reverse transcription pcr)|(rtpcr)|(rt-pcr)|(rt-qpcr)|(rtq-pcr)', re.I)
def pubtools_rtPcrSearch(metaInfo, text):
    """ interface for pubtools """
    #hits = rtPcrRegex.finditer(text)
    match = rtPcrRegex.search(text)
    lines = []
    if match!=None:
        start = match.start()
        end = match.end()
        matchStr = text[start:end]
        #print start, end, matchStr
        fields = [str(start), str(end), matchStr]
        line = "\t".join(fields)
        lines.append(line)
    return "\n".join(lines)

