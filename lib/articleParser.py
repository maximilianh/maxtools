import logging
try:
    from lxml import etree # you need to install this. Debian/Redhat package: python-lxml, see also: codespeak.net/lxml/installation.html
    import lxml
except ImportError:
    import xml.etree.cElementTree as etree
    logging.debug("running with python2.5 cElementTree")

try:
    import lxml.html
    lxmlLoaded = True
except:
    lxmlLoaded = False
    logging.debug("lxml.html not loaded")

import util, namedtuple, html, maxTables, maxXml, articleDownload
import urlparse, time, codecs, os, logging, urllib, urllib2, re, sys, \
    tarfile, zlib, tempfile, shutil, socket, subprocess, random, types
try:
    import unidecode
except:
    pass

PROGNAME = "pubtools"

# the hardcoded data sources (others can be defined in config file)
BUILTIN_DATASETS=["pmcftp", "genetics", "elsevier"]
# 
ELSEVIER_ARTICLE_TAGS = ["converted-article", "article", "simple-article", "book-review", "exam", "book", "chapter", "simple-chapter", "examination", "fb-non-chapter", "glossary", 
"index", "introduction", "bibliography"]

# ====== TEST DATA =====
def PubmedTestDoc():
    return '''
        <PubmedArticle> 
        <MedlineCitation Owner="NLM" Status="MEDLINE">  
        <PMID>20430833</PMID>
        <DateCreated>
            <Year>2010</Year>
            <Month>05</Month>
            <Day>31</Day>
        </DateCreated>
        <DateCompleted>
            <Year>2010</Year>
            <Month>06</Month>
            <Day>22</Day>
        </DateCompleted>
        <Article PubModel="Print-Electronic">
            <Journal>
                <ISSN IssnType="Print">1460-1234</ISSN>
                <ISSN IssnType="Electronic">1460-2156</ISSN>
                <JournalIssue CitedMedium="Internet">
                    <Volume>133</Volume>
                    <Issue>Pt 6</Issue>
                    <PubDate>
                        <Year>2010</Year>
                        <Month>Jun</Month>
                    </PubDate>
                </JournalIssue>
                <Title>Brain : a journal of neurology</Title>
                <ISOAbbreviation>Brain</ISOAbbreviation>
            </Journal>
            <ArticleTitle>Tyrosine hydroxylase deficiency: a treatable disorder of brain catecholamine biosynthesis.</ArticleTitle>
            <Pagination>
                <MedlinePgn>1810-22</MedlinePgn>
            </Pagination>
            <Abstract>
                <AbstractText>An infantile onset, progressive, hypokinetic-rigid syndrome with dystonia (type A), and a complex encephalopathy with neonatal onset (type B). Decreased cerebrospinal fluid concentrations of homovanillic acid and c.698G&gt;A and c.707T&gt;C mutations. Carriership of at least one promotor mutation, however, apparently predicts type A tyrosine hydroxylase deficiency. Most patients with tyrosine hydroxylase deficiency can be successfully treated with l-dopa.</AbstractText>
            </Abstract>
            <Affiliation>Radboud University Nijmegen Medical Centre, Donders Institute for Brain, Cognition and Behaviour, Department of Paediatric Neurology (820 IKNC), PO Box 9101, 6500 HB Nijmegen, The Netherlands. m.willemsen@cukz.umcn.nl</Affiliation>
            <AuthorList CompleteYN="Y">
                <Author ValidYN="Y">
                    <LastName>Willemsen</LastName>
                    <ForeName>Mich&#233;l A</ForeName>
                    <Initials>MA</Initials>
                </Author>
                <Author ValidYN="Y">
                    <LastName>Verbeek</LastName>
                    <ForeName>Marcel M</ForeName>
                    <Initials>MM</Initials>
                </Author> 
            </AuthorList>
            <Language>eng</Language>
            <PublicationTypeList>
                <PublicationType>Journal Article</PublicationType>
                <PublicationType>Research Support, Non-U.S. Gov't</PublicationType>
            </PublicationTypeList>
            <ArticleDate DateType="Electronic">
                <Year>2010</Year>
                <Month>04</Month>
                <Day>29</Day>
            </ArticleDate>
        </Article>
        <MedlineJournalInfo>
            <Country>England</Country>
            <MedlineTA>Brain</MedlineTA>
            <NlmUniqueID>0372537</NlmUniqueID>
            <ISSNLinking>0006-8950</ISSNLinking>
        </MedlineJournalInfo>
        <ChemicalList>
            <Chemical>
                <RegistryNumber>0</RegistryNumber>
                <NameOfSubstance>Catecholamines</NameOfSubstance>
            </Chemical>
        </ChemicalList>
        <CitationSubset>AIM</CitationSubset>
        <CitationSubset>IM</CitationSubset>
        <MeshHeadingList>
            <MeshHeading>
                <DescriptorName MajorTopicYN="N">Age of Onset</DescriptorName>
            </MeshHeading>
            <MeshHeading>
                <DescriptorName MajorTopicYN="Y">Useless Research</DescriptorName>
            </MeshHeading>
        </MeshHeadingList>
        </MedlineCitation>  
        <PubmedData>
           <History>
            <PubMedPubDate PubStatus="aheadofprint">
                <Year>2010</Year>
                <Month>4</Month>
                <Day>29</Day>
            </PubMedPubDate>
            <PubMedPubDate PubStatus="pubmed">
                <Year>2010</Year>
                <Month>5</Month>
                <Day>1</Day>
                <Hour>6</Hour>
                <Minute>0</Minute>
            </PubMedPubDate>
            <PubMedPubDate PubStatus="medline">
                <Year>2010</Year>
                <Month>6</Month>
                <Day>23</Day>
                <Hour>6</Hour>
                <Minute>0</Minute>
            </PubMedPubDate>
        </History>
        <PublicationStatus>ppublish</PublicationStatus>
        <ArticleIdList>
            <ArticleId IdType="pii">awq087</ArticleId>
            <ArticleId IdType="doi">10.1093/brain/awq087</ArticleId>
            <ArticleId IdType="pubmed">20430833</ArticleId>
        </ArticleIdList>
        </PubmedData>
        </PubmedArticle>'''

def which(file, mode=os.F_OK | os.X_OK, path=None):
   """Locate a file in the user's path, or a supplied path. The function
   yields full paths in which the given file matches a file in a directory on
   the path. from bugs.python.org/file15381/
   """
   if not path:
       path = os.environ.get("PATH", os.defpath)

   for dir in path.split(os.pathsep):
       full_path = os.path.join(dir, file)
       if os.path.exists(full_path) and os.access(full_path, mode):
           yield full_path

def downloadFile(url, filename):
    """ download file with curl or wget or python internal to filename """
    if len(list(which("curl")))>0:
        cmdLine = "curl '%s' -s -o '%s'" % (url, filename)
        logging.log(5, "Found curl, running %s" % cmdLine)
    elif len(list(which("wget")))>0:
        cmdLine = "wget '%s' -q -O '%s'" % (url, filename)
        logging.log(5, "Found wget, running %s" % cmdLine)
    else:
        cmdLine = None

    if cmdLine:
        ret = os.system(cmdLine)
        if ret!=0:
            raise FtToolsError("Problem running command %s" % cmdLine)
    else:
        req = urllib2.Request(url)
        opener = urllib2.build_opener(SmartRedirectHandler())   
        req.add_header('User-Agent', 'Mozilla/5.0 (X11; U; Linux i686; en-US; rv:1.7.8) Gecko/20050524 Fedora/1.5 Firefox/1.5')
        fh = opener.open(req)
        open(filename, "wb").write(fh.read())

def setupLogging(logFilename, options):
    """ direct logging to a file PROGNAME.log and also to stdout, depending on options (debug, verbose, nodeId, etc) """
    if options.verbose:
        stdoutLevel=5
        fileLevel=5
    elif options.debug:
        stdoutLevel=logging.DEBUG
        fileLevel=logging.DEBUG
    else:
        stdoutLevel=logging.INFO
        fileLevel=logging.DEBUG

    if options.nodeId!=None:
        logFilename+= "."+str(options.nodeId)

    rootLog = logging.getLogger('')
    rootLog.setLevel(fileLevel)

    logging.root.handlers = []
    logging.basicConfig(level=fileLevel,
                        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                        datefmt='%m-%d %H:%M',
                        filename= logFilename,
                        filemode='w', stream=None)

    #print logging.getLogger('').getEffectiveLevel()

    # define a handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler()
    # set a format which is simpler for console use
    formatter = logging.Formatter('%(levelname)-8s-%(message)s')
    # tell the handler to use this format
    console.setFormatter(formatter)
    console.setLevel(stdoutLevel)
    # add the handler to the root logger
    logging.getLogger('').addHandler(console)

class FtToolsError(Exception):
    pass

class FtExtractError(Exception):
    pass

# --- Meta Info Record and file handling ---
METAINFOFIELDS=[
"numId",  # internal number that identifies this article in the pubtools system
"id", # internal string id of the article, e.g. PMC12343 or doi:123213213/dfsdf or PMID123123
"source",  # the origin of the article, something like "elsevier" or "pubmed" or "medline"
"journal",      # journal or book title
"printIssn",    # ISSN of the print edition of the article
"eIssn",        # optional: ISSN of the electronic edition of the journal/book of the article
"year",         # first year of publication (electronic or print or advanced access)
"articleType", # research-article, review or other
"articleSection",  # the section of the book/journal, e.g. "methods", "chapter 5" or "Comments"
"authors",  # list of authors, often separated by semicolon, but not always
"title",    # title of article
"abstract", # abstract if available
"vol",      # volume
"issue",    # issue
"page",            # first page of article, can be ix, x, or S4
"pmid",            # PubmedID if available
"pmcId",           # Pubmed Central ID
"doi",             # DOI, without leading doi:
"fulltextFileXml", # local filename that contains the full XML
"fulltextFile",    # local filename of ASCII fulltext
"fulltextFilePdf", # local filename of PDF fulltext
"origFile",        # original filename, assigned by publisher
"suppFiles",       # comma-separated list of supplemental filenames
"suppDescs",       # comma-separated list of supplemental descriptions
"origSuppFiles",   # comma-sep. list of original supplementary filenames
"suppUrls",        # comma-seq list of original supp URLs
"fulltextUrl",     # URL to fulltext of article
"downloadDate"     # date of download
]

MetaInfo = namedtuple.namedtuple("MetaInfoRec", METAINFOFIELDS)
emptyMeta = MetaInfo(*len(METAINFOFIELDS)*[""])

def createEmptyMetaDict(pmcId=None, source=None, journal=None, id=None, numId=None, origFile=None):
    metaInfo = emptyMeta._asdict()
    metaInfo["downloadDate"]=time.asctime()
    if pmcId:
        metaInfo["pmcId"]=pmcId
    if origFile:
        metaInfo["origFile"]=origFile
    if source:
        metaInfo["source"]=source
    if journal:
        metaInfo["journal"]=journal
    if id:
        metaInfo["id"]=id
    if numId:
        metaInfo["numId"]=str(numId)
    return metaInfo

def writeHeadersToMetaInfo(path):
    """ create file and add headers to it """
    fh = open(path, "w")
    fh.write("#"+"\t".join(emptyMeta._fields)+"\n")
    fh.close()

def appendMetaInfo(path, metaInfoDict):
    logging.log(5, "appending metaInfo %s to %s" % (str(metaInfoDict), path))
    fh = codecs.open( path, "a", encoding="utf8" )
    if len(metaInfoDict)!=len(METAINFOFIELDS):
        logging.error("column counts between meta file and internal objects don't match")
        metaInfoFields = list(METAINFOFIELDS)
        metaInfoFields.sort()
        logging.error("Internal columns are: %s" % str(metaInfoFields))
        dictKeys = list(metaInfoDict.keys())
        dictKeys.sort()
        logging.error("supplied dict is:     %s" % str(dictKeys))
        sys.exit(1)
    metaInfoTuple = MetaInfo(**metaInfoDict)
    utf8Tuple = []
    for mit in metaInfoTuple:
        if type(mit)==type(1):
            mit = str(mit)
        if mit==None:
            utf8Tuple.append("None")
        elif type(mit)==type(unicode()):
            utf8Tuple.append(mit)
        else:
            utf8Tuple.append(mit.decode("utf8"))
    metaInfoTuple = utf8Tuple
    #metaInfoTuple = [x.decode("utf8")  for x in metaInfoTuple]
    #print metaInfoTuple
    line = "\t".join(metaInfoTuple)
    line = line.replace("\n", " ")
    fh.write(line+"\n")
    fh.close()

# -- convenience
def urlJoin(base, url):
    return urlparse.urljoin(base, url)

def etreeFromHtml(string):
    tree   = lxml.html.document_fromstring(string).xpath("/html/body")[0]
    return tree
    
def etreeFromXml(string):
    """ parses string to etree and removes all namespaces """
    if lxmlLoaded:
        tree   = etree.fromstring(string).getroottree().getroot()
    else:
        tree = etree.fromstring(string)
        #print dir(tree)
        #tree = tree.getroot()
    strip_namespace_inplace(tree)
    return tree

def writeFileObj(fileObj, path):
    logging.debug("Writing to file %s" %path)
    open(path, "w").write(fileObj.read())

def makeAllSubdirs(dirname):
    logging.log(5, "Making all directories for %s" %dirname)
    try:
        os.makedirs(dirname)
    except OSError, ex:
        if ex.errno==17:
            pass
        else:
            raise
    return

def generate2LevelPath(baseDir, filename):
    """ generate a two-level directory path for filename in baseDir, e.g. test.pdf will mkdir $(baseDir)/b1/74/ and return the string $(basedir)/b1/74/test.pdf"""
    hashVal = hex(hash(filename) % 2**16)
    hashVal = hashVal.split("x")[1]
    while len(hashVal)<4:
        hashVal+="0"
    logging.log(5, "Generated directory hash value from filename %s is %s" % (filename, hashVal))
    dir1 = hashVal[:2]
    dir2 = hashVal[2:]
    basedir = os.path.join(baseDir, dir1, dir2)
    makeAllSubdirs(basedir)
    return os.path.join(basedir, filename)

# --- retrieve Indexes

def queryEntrez(dbName, query, onlyCounts=False, minYear=None, maxYear=None, retStart=0):
    """ retrieve pmids for query, returns a tuple (no of results, list of ids) optionally returns only number of results """
    if onlyCounts:
        retmax=0
    else:
        retmax=9999999

    if minYear!=None and maxYear!=None:
        addString = "&mindate=%s&maxdate=%s" % (minYear, maxYear)
    else:
        addString = ""

    query = urllib.quote(query)
    url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=%s&tool=retrPubmed&email=maximilianh@gmail.com&term=%s&retstart=%d&retmax=%d%s' % (dbName, query, retStart, retmax, addString)
    req = urllib2.Request(url)
    req.add_header('User-Agent', 'Mozilla/5.0 (X11; U; Linux i686; en-US; rv:1.7.8) Gecko/20050524 Fedora/1.5 Firefox/1.5')
    html = urllib2.urlopen(req)
    logging.info("Retrieving "+url+"\n")
    
    pmids = []
    count = None
    for line in html:
        if line.find("<Id>")!=-1:
            pmid = line.strip().replace("<Id>","").replace("</Id>", "")
            pmids.append(pmid)
        if line.find("<Count>")!=-1:
            if count==None: # only use the first "count" line
                fs = re.split("<.{0,1}Count>", line)
                count = int(fs[1])
                if count == "-1":
                    logging.info("No ids found")
                    return []

    logging.info("Query matches %s records. Retrieved %s records from Pubmed. \n" % (str(count), str(len(pmids))))

    if int(count)!=len(pmids):
        logging.warn("Warning: not all data was retrieved. You need to set the retstart parameter and run this script again\n")

    return pmids

# ---- downloaders -----------
def downloadOneGeneticsArticle(pmid, outDir, metaInfoFilename, minNumId):
        pmid = pmid.strip()
        logging.info("Getting urls for PMID"+pmid)
        artId = "genetics-PMID"+str(pmid)

        metaInfo = createEmptyMetaDict()
        metaInfo["pmid"]=pmid
        metaInfo["source"]="genetics"
        metaInfo["journal"]="Genetics"
        metaInfo["id"]=artId

        urlDict = pmidToGeneticsUrls(pmid)
        #print urlDict
        filestem = generate2LevelPath(outDir, artId)
        logging.log(5,"Got URLs from genetics: %s " % str(urlDict))

        if "pdf" in urlDict:
            filename = filestem+".pdf"
            util.httpDownload(urlDict["pdf"], filename, verbose=True)
            metaInfo["fulltextFilePdf"] = util.relpath(filename, outDir)
            metaInfo["fulltextUrl"] = urlDict["pdf"].replace(".pdf","")

        if "full" in urlDict:
            filename = filestem+".html"
            util.httpDownload(urlDict["full"], filename, verbose=True)
            metaInfo["fulltextFile"] = util.relpath(filename, outDir)
            metaInfo["fulltextUrl"] = urlDict["full"]

        if "suppData" in urlDict:
            i=1
            suppnames = []
            suppurls = []
            suppFiles = []
            for desc, suppUrl in urlDict["suppData"]:
                filext = suppUrl.split(".")[-1]
                filename = filestem+".supp%d.%s" % (i, filext)
                util.httpDownload(suppUrl, filename, verbose=True)
                suppnames.append(filename)
                suppurls.append(suppUrl)
                suppFiles.append(util.relpath(filename, outDir))
                i+=1
            metaInfo["suppUrls"]=",".join(suppurls)
            metaInfo["suppDescs"]=",".join(suppnames)
            metaInfo["origSuppFiles"]=",".join(suppFiles)
            metaInfo["numId"]=(minNumId+1)
        appendMetaInfo(metaInfoFilename, metaInfo) 

def downloadGenetics(pmids, outDir, metaInfoFilename, sleepSeconds, minNumId, testId):
    pass

def extractPmcFile(pmcBaseUrl, tarfname, tempExtractDir):
    """ extract all files from pmc tarfile, returns tuple of (filenames, dirnames)"""

    # open and extract tarfile into temp dir
    if not os.path.exists(tarfname):
        logging.warn("Could not find file %s, internal error?\n" % tarfname)
        return [], []
    logging.debug("Opening tarfile %s" % tarfname)
    try:
        tarObject = tarfile.open(tarfname)
    except:
        raise FtExtractError("error while while opening tar file")

    if tarObject==None:
        raise FtToolsError("Could not find %s: index file not in sync with pmc ftp folder" % (tarfname))
        return [], []

    logging.debug("Extracting tarfile %s to %s" % (tarfname, tempExtractDir))
    # tarObject.extractall(tempExtractDir) # only in py2.5

    try:
        util.extractTar(tarObject, tempExtractDir)
    except:
        raise FtExtractError("error while extracting from tar file")

    #cmdLine = "tar xfz %s -C %s" % (tarfname, tempExtractDir)
    #ret = os.system(cmdLine)
    #if ret!=0:
    #    raise FtExtractError("Error while untarring %s" % (tarfname))
        
    filenames = []
    dirnames = []
    for m in tarObject.getmembers():
        if m.isfile():
            filenames.append(os.path.join(tempExtractDir, m.name))
        if m.isdir():
            dirnames.append(os.path.join(tempExtractDir, m.name))

    return filenames, dirnames

def processTempPmcArticleFolder(filenames, metaInfoFilename, pmcId, repoDir, numId, origFile):
    """ extract metaInfo from pmc article directory and copy files into repository """

    #filenames = util.findSubdirFiles(path, None)
    nxmlfiles = [f for f in filenames if f.endswith(".nxml")]
    if len(nxmlfiles)!=1:
        raise FtToolsError("Illegal number of nxml files in files %s, pmcId %s" % (str(filenames), pmcId))
    metaInfo = createEmptyMetaDict(source = "pmc", id = pmcId, pmcId = pmcId.replace("PMC", ""), numId=numId, origFile=origFile)
    nxmlName = nxmlfiles[0]
    xmlData, metaInfo = parsePmcXml(nxmlName, metaInfo)
    if metaInfo["articleType"]=="noContent":
        logging.warn("Article XML includes no <article>-tag, skipping further import")
        appendMetaInfo(metaInfoFilename, metaInfo)
        return 

    pmcId    = metaInfo["id"]

    filestem = generate2LevelPath(repoDir, pmcId)

    localNxmlName = filestem+".nxml"
    logging.debug("Main text : Writing xml data for article with pmcId %s to %s" % (pmcId, localNxmlName))
    open(localNxmlName,"w").write(xmlData)
    metaInfo["fulltextFileXml"]=util.relpath(localNxmlName, repoDir)

    if nxmlName!=None:
        nxmlBase = os.path.splitext(nxmlName)[0]
        for altExt, metaTag in [('.pdf', 'fulltextFilePdf'), ('.txt', 'fulltextFile')]:
            filename = nxmlBase+altExt
            logging.log(5, "checking if file %s exists" % (filename))
            if filename in filenames:
                extension = os.path.splitext(filename)[1]
                newName  = filestem+extension.lower()
                logging.debug("Main Text : Copying %s to %s" % (filename, newName))
                shutil.copy(filename, newName)
                metaInfo[metaTag] = util.relpath(newName, repoDir)
    else:
        logging.warn("filenames %s do not contain an nxml file!" % str(filenames))

    # copy supp files
    baseDir = os.path.dirname(nxmlName)
    i = 1
    suppFiles = metaInfo["origSuppFiles"]
    newSuppNames = []
    logging.log(5, "suppFiles entry is %s" % suppFiles)
    if suppFiles!=None and suppFiles!="":
        for suppFile in suppFiles.split(","):
            origName = os.path.join(baseDir, suppFile) 
            ext = os.path.splitext(suppFile)[1]
            newName = filestem+".S"+str(i)+ext.lower()
            logging.log(5, "Supp Data : Copying %s to %s" % (origName, newName))
            shutil.copy(origName, newName)
            newSuppNames.append(util.relpath(newName, repoDir))
            i+=1
    metaInfo["suppFiles"]=",".join(newSuppNames)

    metaInfo["numId"] = numId
    appendMetaInfo(metaInfoFilename, metaInfo)

def deleteTmpFilesDirs(dirnames):
    if len(dirnames)==0:
        return

    dirnames.sort(key=len)
    shortestDirname = dirnames[0]
    if shortestDirname!=None:
        logging.debug("Deleting with all subdirectories: %s" % shortestDirname)
        shutil.rmtree(shortestDirname)

def downloadToCache(localTarName, pmcBaseUrl, tarfname, localCacheDir):
    tarUrl = urlJoin(pmcBaseUrl, tarfname)
    localDir = os.path.dirname(localTarName)
    if not os.path.isdir(localDir):
        logging.debug("creating directory %s" % localDir)
        os.makedirs(localDir)
    else:
        logging.log(5, "directory %s already exists" % localDir)

    logging.info("Downloading %s to %s" % (tarUrl, localTarName))
    #tarData = urllib2.urlopen(tarUrl).read()
    #locTarFile = open(localTarName, "w")
    #locTarFile.write(tarData)
    #locTarFile.close()
    downloadFile(tarUrl, localTarName)
    logging.debug("Wrote data to %s" % localTarName)

def downloadOnePmc(pmcId, pmcBaseUrl, tarfname, metaInfoFilename, tempExtractDir, localCacheDir, repoDir, numId, doDownload, testMode=False):

    filenames = []
    toDelDirnames = []

    try:
        localTarName = os.path.join(localCacheDir, tarfname)
        if (not os.path.isfile(localTarName)) or os.path.getsize(localTarName)==0:
            if pmcBaseUrl!=None and doDownload:
                logging.debug("File %s not found in localCacheDir, need to download it" % localTarName)
                downloadToCache(localTarName, pmcBaseUrl, tarfname, localCacheDir)
            elif doDownload:
                raise FtToolsError("File %s not found, but pmcBaseUrl not set, skipping file" % localTarName)
            else:
                logging.info("File %s not found, but doDownload not set, skipping this file" % localTarName)
                return

        try:
            filenames, toDelDirnames = extractPmcFile(pmcBaseUrl, localTarName, tempExtractDir)
        # we handle some exceptions by re-downloading + retrying
        except (IOError, FtExtractError):
            logging.warn("Exception (%s) on extraction, trying to re-download and continue" %sys.exc_info()[1])
            if pmcBaseUrl==None:
                raise FtToolsError("pmcBaseUrl not set, skipping file" )
            else:
                downloadToCache(localTarName, pmcBaseUrl, tarfname, localCacheDir)
                filenames, toDelDirnames = extractPmcFile(pmcBaseUrl, localTarName, tempExtractDir)
            
        try:
            processTempPmcArticleFolder(filenames, metaInfoFilename, pmcId, repoDir, numId, tarfname)
        except lxml.etree.XMLSyntaxError, SyntaxError:
            logging.warn("Exception (%s) on XML parsing, skipping this article" %sys.exc_info()[1])
        deleteTmpFilesDirs(toDelDirnames)

    # these are exceptions that completely skip the article
    except Exception, ex:
        if testMode:
            logging.debug("Exception %s occured, but testing mode, keeping all files"%sys.exc_info()[1])
            raise
        else:
            logging.debug("Exception %s occured, completely zapping temp dir, skipping this article"%sys.exc_info()[1])
            shutil.rmtree(tempExtractDir, ignore_errors=True)
            os.makedirs(tempExtractDir)

def runFunctionOnTextfile(indexFile, processLine, minNumId, nodeCount, nodeId, testId, skipIds):
    """ run a function on lines of a textfile. Use only each nodeCount line, starting from nodeId.
        supply the function with testId (skip all but this ID) and skipIds (skip these IDs)
    """

    # make sure that numIds do not overlap between two nodes
    if nodeId==None:
        idStepSize=1
        numId = minNumId+1
    else:
        idStepSize=nodeCount
        numId = minNumId+1+int(nodeId)

    lineNo = 0
    for line in indexFile:
        line = line.strip()
        lineNo+=1
        if nodeCount!=None and lineNo%nodeCount!=nodeId:
            logging.log(5, "Skipping %s as it is not for this cluster node" % line)
            continue
        else:
            if lineNo % 10000 == 0 and not testId:
                logging.info("Processed %d entries" % lineNo)
        if testId!=None:
            processLine(line, testId, numId,skipIds) # in testing mode, we do not want to catch exceptions
        else:
            try:
                processLine(line, testId, numId,skipIds)
            except IOError, ex:
                if ex.errno==28:
                    logging.error("Disk seems to be full, stopping processing completely")
                    raise
                else:
                    logging.warn("IO Exception %s: %s on line %s" % (sys.exc_info()[0], sys.exc_info()[1], line))
            except (FtExtractError, FtToolsError):
                logging.warn("Exception %s: %s, on line %s" % (sys.exc_info()[0], sys.exc_info()[1], line))
            numId+=idStepSize

def readOldMetaFile(metaInfoFilename, minNumId, idField="id"):
    """ read the fields "numId" and "id" from old metaInfoFile and return 
        the highest numId and a set of all string ids found """
    doneIds = []
    minNumId = int(minNumId)
    if os.path.isfile(metaInfoFilename):
        logging.debug("Reading processed doneIds from metaInfo file %s" % (metaInfoFilename))
        numIdsAndIds = list(maxTables.TableParser(open(metaInfoFilename)).columns(["numId", idField]))
        doneIds = set([y for (x,y) in numIdsAndIds])
        logging.info("Read %d doneIds from metaInfo file (=should be already in the repository)" % len(doneIds))
        numIds = set([int(x) for (x,y) in numIdsAndIds])
        if len(numIds)!=0:
            minNumId = max(numIds)
    return doneIds, minNumId

def download_Pmc(indexFilename, metaInfoFilename, pmcBaseUrl, tempExtractDir, localCacheDir, repoDir, minNumId, nodeCount, nodeId, doDownload, testId=None):
    """ download and parse pmc from ftp server to localCacheDir, skip files
    that are already in localCacheDir, skip files that are already in metaInfo
    """
    # little convenience function
    def processLine(line, testId, numId, skipIds):
        tarfname,journal,pmcId = line.strip().split("\t")
        if testId:
            if pmcId!=testId:
                return
        if pmcId in skipIds and not testId:
            logging.debug("Skipping %s, already in metainfo file" % (pmcId))
            return
        logging.debug("====Processing %s, tarfile %s====" % (pmcId, tarfname))
        testMode = (testId!=None)
        downloadOnePmc(pmcId, pmcBaseUrl, tarfname, metaInfoFilename, tempExtractDir, localCacheDir, repoDir, numId, doDownload, testMode=testMode)

    # FUNCTION START
    logging.debug("Cleaning Temp directory %s" % tempExtractDir)
    shutil.rmtree(tempExtractDir, ignore_errors=True)
    os.makedirs(tempExtractDir)

    fileListPath = os.path.join(localCacheDir, "file_list.txt")
    if not os.path.isfile(fileListPath):
        logging.warn('Sure that "localCacheDir=" is correct in pubtools.conf? Could not find %s' % fileListPath) 

    skipIds, minNumId = readOldMetaFile(metaInfoFilename, minNumId)
    logging.debug("starting numerical ID is %d" % minNumId)

    indexFile = open(indexFilename)
    indexFile.readline() # skip header of pmc file

    runFunctionOnTextfile(indexFile, processLine, minNumId, nodeCount, nodeId, testId, skipIds)


# --- parsers ---

def metaInfoFromTree(metaInfo, field, etree, xpath, convertToAscii=False):
    """ helper function for parsePmcXml, get text given xpath expression, put into metaInfo[field] """ 
    elements = etree.xpath(xpath)

    if len(elements)==0:
        return
    if len(elements)>1:
        #for e in elements:
            #print e.attrib, e.text
        logging.warn("xpath expression %s lead to more than one value, only first one is used" % xpath)

    element = elements[0]

    if convertToAscii==True:
        metaInfo[field] = treeToAsciiText(element)
    else:
        metaInfo[field] = element.text

    if metaInfo[field] != None:
        metaInfo[field] = metaInfo[field].replace("\n", " ")
        metaInfo[field] = metaInfo[field].replace("\t", " ")


def extractNames(contribNodeList):
    """ helper function for parsePmcXml, get the author names from the contrib
    section of a pmc xml file""" 
    surnames = []
    firstnames = []
    for contrib in contribNodeList:
        contribType = contrib.attrib.get("contrib-type",None)
        #print etree.tostring(contrib, pretty_print=True)
        if contribType=="author":
            surnameEls = contrib.xpath("name/surname")
            if len(surnameEls)!=0:
                surname = surnameEls[0].text
                if surname==None:
                    surname=""
                surnames.append(surname)
            firstnameEls = contrib.xpath("name/given-names")
            if len(firstnameEls)!=0:
                firstname = firstnameEls[0].text
                if firstname==None:
                    firstname=""
                firstnames.append(firstname)
    names = []
    for i in range(0, min(len(surnames), len(firstnames))):
        names.append(surnames[i]+", "+firstnames[i])
    return "; ".join(names)

def minimumYear(articleMeta):
    """ helper for parsePmcXml: extracts the minimum year of all pubdate entries """
    pubDates = articleMeta.xpath("pub-date")
    years = []

    if pubDates!=None:
        for date in pubDates:
            year = date.xpath("year")[0]
            if year==None:
                continue
            yearString = year.text
            try:
                years.append(int(yearString))
            except:
                logging.debug("exception when casting the year to int")
                pass
        
    if len(years)==0:
        logging.debug("no valid years found")
        return "noValidYear"

    minYear = str(min(years))
    logging.debug("determined minimum year as %s" % minYear)
    return minYear

def treeToAsciiText(tree, _addtail=False):
    return "".join(recursiveToAscii(tree))

def recursiveToAscii(tree, _addtail=True):
    """ xml -> ascii tags: convert all text associated with all tags to a
    space-sep. ASCII text string in utf8 
    copied from http://code.activestate.com/recipes/498286/ 
    Returns a list of text strings contained within an element and its sub-elements.
    Helpful for extracting text from prose-oriented XML (such as XHTML or DocBook).
    """
    result = []
    if tree.text is not None:
        result.append(tree.text)
    for elem in tree:
        result.extend(recursiveToAscii(elem,True))
    if _addtail and tree.tail is not None:
        result.append(tree.tail)
    return result

def parsePmcXml(nxmlName, metaInfo):
    """ return the meta information contained in an pmc xml file as a metaInfo-dictionary """

    logging.debug("Parsing pmc nxml file %s" % nxmlName)
    #tree = etreeFromHtml(xml) # got inserted automatically by html parser
    #xml  = codecs.open(nxmlName, encoding="utf8").read()
    xml  = open(nxmlName).read()
    articleTree = etreeFromXml(xml) # got inserted automatically by html parser
    #print etree.tostring(tree, method="text", encoding="utf8")

    if articleTree==None or len(articleTree)==0 or articleTree.tag!="article":
        if articleTree.tag!="article":
            logging.info("Not an article: %s tag found" % articleTree.tag)
        else:
            logging.warn("No article tag found, empty file?")

        metaInfo["articleType"]="noContent"
        return xml, metaInfo

    metaInfo["articleType"]=articleTree.xpath("@article-type")[0]
    journalTree = articleTree.xpath("front/journal-meta")[0]
    metaInfoFromTree(metaInfo, "journal", journalTree, "journal-title")
    if metaInfo["journal"]=="":
        metaInfoFromTree(metaInfo, "journal", journalTree, "journal-title-group/journal-title")
    metaInfoFromTree(metaInfo, "printIssn", journalTree, "issn[@pub-type='ppub']")
    metaInfoFromTree(metaInfo, "eIssn", journalTree, "issn[@pub-type='epub']")

    articleTree = articleTree.xpath("front/article-meta")[0]
    metaInfoFromTree(metaInfo, "pmcId", articleTree, "article-id[@pub-id-type='pmc']")
    metaInfo["id"]="PMC"+metaInfo["pmcId"]
    metaInfoFromTree(metaInfo, "pmid", articleTree, "article-id[@pub-id-type='pmid']")
    metaInfoFromTree(metaInfo, "title", articleTree, "title-group/article-title", convertToAscii=True)
    if metaInfo["title"]==None or metaInfo["title"]=="":
        metaInfoFromTree(metaInfo, "title", articleTree, "title-group/article-title", convertToAscii=True)

    minYear = minimumYear(articleTree)
    metaInfo["year"]=minYear

    contribs = articleTree.xpath("contrib-group/contrib")
    nameStr = extractNames(contribs)
    metaInfo["authors"]=nameStr

    metaInfoFromTree(metaInfo, "vol", articleTree, "volume")
    metaInfoFromTree(metaInfo, "issue", articleTree, "issue")
    metaInfoFromTree(metaInfo, "page", articleTree, "fpage")
    metaInfoFromTree(metaInfo, "abstract", articleTree, "abstract", convertToAscii=True)
    #print etree.tostring(tree, pretty_print=True)

    # supp data files
    suppFiles = []
    suppDescs = []
    for suppData in articleTree.xpath("//supplementary-material"):
        shortDesc = suppData.attrib.get("id","")
        longDescList = suppData.xpath("caption/title")
        longDesc = ""
        if len(longDescList)>0:
            longDesc = longDescList[0].text
        if longDesc==None:
            longDesc=""
        if len(longDesc)>len(shortDesc):
            desc = longDesc
        else:
            desc = shortDesc
        desc = desc.replace(",",";") # we use , as a field separator
        medias = suppData.xpath("media")
        for media in medias:
            fname = media.attrib.get("href","")
            suppFiles.append(fname)
            suppDescs.append(desc)

    metaInfo["origSuppFiles"]=",".join(suppFiles)
    metaInfo["suppDescs"]=",".join(suppDescs)

    # remove all tabs and returns from metaInfo
    newMetaInfo = {}
    for key, value in metaInfo.iteritems():
        if value!=None:
            value = value.replace("\n", " ")
            value = value.replace("\t", " ")
        newMetaInfo[key]=value
    metaInfo=newMetaInfo

    logging.log(5, "Complete Meta info:"+str(metaInfo))
    return xml, metaInfo
    
def pmidToGeneticsUrls(pmid):
    """ go to the genetics website and get the urls to the genetics article as a dictionary """

    pmid = pmid.strip()
    url = "http://www.genetics.org/cgi/pmidlookup?view=short&pmid=%s" % pmid
    logging.log(5, "Getting %s" % url)
    htmlFile = util.httpGet(url)
    html = htmlFile.read()
    baseUrl = htmlFile.url

    urls = {}
    if "reprint" in baseUrl:
        logging.log(5, "looks like direct link to pdf, no parsing of html")
        urls["pdf"]=baseUrl+".pdf"
    else:
        linkDict = {} 
        tree = etreeFromHtml(html)
        links = tree.xpath("div/div/div/div/ul/li/a")
        for l in links:
            #print (l.tag, l.text_content().strip(), l.tail, l.attrib["href"])
            desc = l.text_content().strip("\t\n")
            linkDict[desc] = l.attrib["href"]

        logging.log(5, "links on page: %s" % str(linkDict))

        urls = {}
        if "Full Text" in linkDict:
            urls["full"] = urlJoin(baseUrl, linkDict.get("Full Text", ""))
        if "Full Text (PDF)" in linkDict:
            urls["pdf"] = linkDict.get("Full Text (PDF)", None)
        if "Full Text (Rapid PDF)" in linkDict:
            urls["pdf"] = linkDict.get("Full Text (Rapid PDF)")

        if "pdf" in urls:
            urls["pdf"] = urlJoin(baseUrl, urls["pdf"]+".pdf")

        logging.log(5, "links on page that were recognized: %s" % str(urls))
        suppUrl = linkDict.get("Data Supplement", None)
    
        if suppUrl!=None:
            suppFile = util.httpGet(urlJoin(baseUrl,suppUrl))
            html = suppFile.read()
            tree = etreeFromHtml(html)
            suppUrls = tree.xpath("div/div/div/div/font/ul/li/a")

            suppUrlList = []
            for l in suppUrls:
                #print (l.tag, l.text_content().strip(), l.tail, )
                suppUrlList.append( (l.text_content().strip(), urlJoin(baseUrl, l.attrib["href"])) )

            urls["suppData"]=suppUrlList
            #print urls
    return urls

def readMetaInfos(filename, nodeCount=None, nodeId=None, skipIds=set()):
    """ generate metaInfo objects from metaFile """
    skipIds = set(skipIds)
    tp = maxTables.TableParser(open(filename))

    i = 0
    for l in tp.lines():
        i+=1
        if l.id in skipIds:
            logging.log(5, "skipping %s, already in skipIds" % l.id)
            continue
        if nodeCount!=None:
            if i%nodeCount==nodeId:
                yield l
        else:
            yield l

def writeMetaInfoDicts(filename, metaInfoDicts):
    """ convert a list of dictionares to MetaInfo-objects and append them to a file """
    writeHeadersToMetaInfo(filename)
    ofh = open(filename, "a")

    for metaInfoDict in metaInfoDicts.values():
        metaInfoTuple = MetaInfo(**metaInfoDict)
        ofh.write("\t".join(metaInfoTuple)+"\n")
    ofh.close()

def toAscii(inFile, outFile, converters, publisher):
    """ convert inFile to Ascii and write to outfile, return False if not
    successful, using converters dictionary.  
    key = extension, value = command to execute, with vars $in and $out """
    baseName, ext = os.path.splitext(inFile)
    ext = ext.strip(".").lower()
    if ext not in converters:
        logging.debug("Could not convert file %s, no converter for extension %s" % (inFile, ext))
        return False
    else:
        cmdLine = converters[ext]
        if os.path.exists(outFile):
            logging.debug("output file %s exists, skipping" % outFile)
            return True
        outDir = os.path.dirname(outFile)
        if not os.path.isdir(outDir):
            try:
                os.makedirs(outDir)
            except OSError, ex:
                if ex.errno==17:
                    print ex
                    pass
                else:
                    raise

        if cmdLine=="COPY":
            logging.debug("Copying %s to %s" % (inFile, outFile))
            try:
                shutil.copy(inFile, outFile)
            except IOError:
                logging.debug("Exception occured while copying")
                return False
            return True
        elif cmdLine=="XMLTEXT" or cmdLine=="NXMLTEXT":
            logging.debug("stripping XML tags from %s and writing to %s" % (inFile, outFile))
            #inData = codecs.open(inFile, encoding="utf8").read()
            inData = open(inFile).read()
            if cmdLine=="NXMLTEXT":
                asciiData = stripXmlTags(inData, isNxmlFormat=True)
            if publisher=="elsevier":
                asciiData = stripXmlTags(inData, isElsevier=True)
            else:
                asciiData = stripXmlTags(inData)
            if asciiData!=None:
                codecs.open(outFile, "w", "utf8").write(asciiData)
            else:
                return False
        else:
            cmdLine = cmdLine.replace("$in", inFile)
            cmdLine = cmdLine.replace("$out", outFile)
            logging.debug("running "+cmdLine)
            ret = os.system(cmdLine)
            if ret==2:
                logging.info("stopped on errno 2: looks like you pressed ctrl-c")
                sys.exit(2)
            if ret!=0:
                logging.warn("error %d occured while executing %s" % (ret, cmdLine))
                return False
            return True

def metaInfoFilenameGenerator(metaInfo):
    """
    return files of metaInfo entry as tuples in format (fieldname, filename)

    >>> mid = createEmptyMetaDict() 
    >>> mid["fulltextFile"]="test.pdf"
    >>> mid["suppFiles"]="test.s1.pdf,test.s2.pdf"
    >>> mi = MetaInfo(**mid)
    >>> print list(metaInfoFilenameGenerator(mi))
    [('fulltextFile', 'test.pdf'), ('suppFiles', 'test.s1.pdf'), ('suppFiles', 'test.s2.pdf')]
    """
    for field in ['fulltextFile', 'fulltextFilePdf', 'fulltextFileXml']:
        data = metaInfo._asdict()[field]
        if data!="":
            yield field, data

    i = 0
    suppFilenames=metaInfo.suppFiles
    if suppFilenames!="":
        suppFilenames = suppFilenames.split(",")
        for suppFile in suppFilenames:
            yield "suppFiles", suppFile
            i+=1

def commaSepRemove(string, index):
    """ remove value at index from comma sep string """
    if "," in string:
        valList = string.split(",")
        #print valList, index
        if index<len(valList):
            del valList[index]
        else:
            logging.error("Inconsistency between suppFiles and another comma-spec field found, index is %d" % index)
        #index = valList.index(value)
        #valList.remove(value)
        return ",".join(valList)
    else:
        if index==0:
            return ""
        else:
            logging.warn("Cannot delete index %d from comma sep string '%s'" % (index, string))
            return string

def removeMetaInfoFilename(metaInfoDict, fieldname, filename):
    """
    remove filename from metaInfoDict's entry of 'fieldname'. 

    >>> mid = {}
    >>> mid["fulltextFile"]="test.pdf"
    >>> mid["suppFiles"]="test.s1.pdf,test.s2.pdf"
    >>> mid["suppDescs"]="S1,S2"
    >>> mid["origSuppFiles"]="tests1.pdf,tests2.pdf"
    >>> mid["suppUrls"]="http://www.tests1.com.pdf,http://www.tests2.com.pdf"
    >>> removeMetaInfoFilename(mid, "fulltextFile", "test.pdf")
    {'suppFiles': 'test.s1.pdf,test.s2.pdf', 'origSuppFiles': 'tests1.pdf,tests2.pdf', 'fulltextFile': '', 'suppUrls': 'http://www.tests1.com.pdf,http://www.tests2.com.pdf', 'suppDescs': 'S1,S2'}
    >>> removeMetaInfoFilename(mid, "suppFiles", "test.s1.pdf")
    {'suppFiles': 'test.s2.pdf', 'origSuppFiles': 'tests2.pdf', 'fulltextFile': '', 'suppUrls': 'http://www.tests2.com.pdf', 'suppDescs': 'S2'}
    """
    if fieldname=="suppFiles":
        index                         = metaInfoDict["suppFiles"].split(",").index(filename)
        metaInfoDict["suppFiles"]     = commaSepRemove(metaInfoDict["suppFiles"], index)
        metaInfoDict["suppDescs"]     = commaSepRemove(metaInfoDict["suppDescs"], index)
        metaInfoDict["origSuppFiles"] = commaSepRemove(metaInfoDict["origSuppFiles"], index)
        metaInfoDict["suppUrls"]      = commaSepRemove(metaInfoDict["suppUrls"], index)
    else:
        metaInfoDict[fieldname]=""
    return metaInfoDict


def generateFilenames(baseDir, publisher):
    """ open metainfo file, parse it and generate many tuples (metaInfo, filenames) """
    for metaInfo in readMetaInfos(os.path.join(baseDir, publisher+".metaInfo")):
        metaDict = metaInfo._asdict()
        fulltextFields = "fulltextFile,fulltextFileXml,fulltextFilePdf".split(",")
        for field in fulltextFields:
            filename = metaDict[field]
            if filename=="":
                continue
            yield metaInfo, filename

        suppFiles = metaInfo.suppFiles
        if suppFiles!="":
            suppFiles = suppFiles.split(",")
            for suppFile in suppFiles:
                yield metaInfo, suppFile
                

def generateText(asciiDir, publisher, nodeCount, nodeId, textType, testId=None):
    """ return tuples of (metaInfo, isSupplemental, filename, textString)
    If textType="onlybest" then yield only one file per article (try, in this order: ascii, pdf, xml)
    If textType="bestsupp": like best, but also all supp files
    If textType="all" then yield all files per article (xml, ascii,pdf + supp files)
    """
    for metaInfo in readMetaInfos(os.path.join(asciiDir, publisher+".metaInfo"), nodeCount, nodeId):
        if testId!=None and metaInfo.id!=testId:
            continue
        fulltextFiles = []
        if metaInfo.fulltextFile!="":
            fulltextFiles.append(metaInfo.fulltextFile)
        if metaInfo.fulltextFilePdf!="":
            fulltextFiles.append(metaInfo.fulltextFilePdf)
        if metaInfo.fulltextFileXml!="":
            fulltextFiles.append(metaInfo.fulltextFileXml)

        if textType=="onlybest" or textType=="best" or textType=="bestsupp":
            if len(fulltextFiles)==0:
                logging.info("No fulltext file found, skipping article."+metaInfo.id)
                continue
            else:
                fulltextFiles = [fulltextFiles.pop(0)]
                logging.debug("processing only best file: %s" % fulltextFiles[0])
            
        for fulltextFilename in fulltextFiles:
            relFulltextFilename = fulltextFilename
            fulltextFilename = os.path.join(asciiDir, fulltextFilename)
            if os.path.isfile(fulltextFilename):
                fulltextString = open(fulltextFilename).read()
            else:
                logging.debug("fulltext file referenced by metaInfo-file does not exist: %s, skipping this file" % fulltextFilename)
                continue
            yield metaInfo, False, relFulltextFilename, fulltextString

        if metaInfo.suppFiles!="" and (textType=="all" or textType=="bestsupp"):
            for relSuppFilename in metaInfo.suppFiles.split(","):
                suppFilename = os.path.join(asciiDir, relSuppFilename)
                if os.path.isfile(suppFilename):
                    suppString = open(suppFilename).read()
                else:
                    logging.log(5, "supplementary file referenced by metaInfo-file does not exist: %s, skipping this file" % suppFilename)
                    suppString=""
                yield metaInfo, True, relSuppFilename, suppString


def mergePubmedIntoMetaInfo(pubmedData, metaInfoDict):
    for field in set(metaInfoDict).intersection(pubmedData):
        string = pubmedData[field] 
        if string!="":
            metaInfoDict[field] = string 
    return metaInfoDict

def pullOutConfig(datasetName, datasetConfig):
    """ return downloader, downloadParams, numIdStart from 
        datasetConfig, given a dataset name  
    """
    configDict = {}
    for name, confString in datasetConfig:
        fs = confString.split(";")
        configDict[name] = (fs[0], fs[1], fs[2]) # downloader, downloadParams/query, numIdStart

    if datasetName not in configDict:
        logging.error("Dataset %s is not defined in config file (defined are: %s and %s)" % (datasetName, ",".join(configDict.keys()), ", ".join(BUILTIN_DATASETS)))
        sys.exit(1)

    datasetName = datasetName.lower()
    downloader, downloadParams, numIdStart = configDict[datasetName]
    if downloader not in ["pubmed", "pubmedcentral"]:
        logging.error("Error in config file: Downloader for dataset %s has to be pubmed or pubmedcentral" % downloader)
        sys.exit(1)

    return downloader, downloadParams, numIdStart

def getIndex_Dataset(publisher, datasetConfig, downloadDir, tool, email):
    """ download index. 
    The index file is actually a .metaInfo file downloaded from Pubmed.
    
    datasetConfig: datasetName = downloader;query/params
    """

    # get config
    downloader, query, numId = pullOutConfig(publisher, datasetConfig)
    numId = int(numId)

    indexFilename = os.path.join(downloadDir, publisher+".index")

    logging.info("Downloading PMIDs that match query %s via EUtils" % query)
    pmids = ncbiESearch(query, tool=PROGNAME, email=email)
    random.shuffle(pmids)
    logging.info("Got %d PMIDs" % (len(pmids)))
    #open(indexFilename, "w").write("\n".join(pmids))
    #logging.info("Written %d PMIDs to file %s" % (len(pmids), indexFilename))

    #downloadType, query, numIdStart = pullOutConfig(publisher, datasetConfig)
    #indexFilename = os.path.join(downloadDir, publisher+".index")

    # find how many pmids' metainformation needs to be downloaded
    #indexPmids    = set([x.strip() for x in open(indexFilename).readlines()])
    #logging.info("Read %d ids from index file" % len(indexPmids))
    #metaInfoPmids = set(maxTables.TableParser(open(metaInfoFilename)).column("pmid"))
    #logging.info("Read %d ids from metaInfo file %s" % (len(metaInfoPmids), metaInfoFilename))
    #toDownloadIds = indexPmids.difference(metaInfoPmids)
    #logging.info("There are %d ids left to download" % len(toDownloadIds))

    # find first valid numId
    #oldNumIds     = set(maxTables.TableParser(open(metaInfoFilename)).column("numId", dataType=types.IntType))
    #if len(oldNumIds)==0:
        #numId = int(numIdStart)
    #else:
        #numId = max(oldNumIds)+1

    writeHeadersToMetaInfo(indexFilename) 
    logging.info("Getting meta data from pubmed")
    #if len(toDownloadIds)>0:
    lines = 0
    for pubmedData in ncbiEFetchGenerator(pmids):
        logging.log(5, "Raw parsed data: %s" % str(pubmedData))
        metaInfoDict = createEmptyMetaDict()
        metaInfoDict = mergePubmedIntoMetaInfo(pubmedData, metaInfoDict)
        metaInfoDict["numId"] = numId
        metaInfoDict["source"] = "pubmedCrawler"
        appendMetaInfo(indexFilename, metaInfoDict)
        lines+=1
        numId+=1
    logging.info("Wrote %d lines to %s" % (lines, indexFilename))
    #else:
        #logging.info("No need to download any data from Pubmed")

def download_Dataset(publisher, datasetConfig, downloadDir, fulltextWriter):
    """ download pdf files via outlinks from pubmed """
    logging.info("Downloading PDFs")
    downloader = articleDownload.FulltextDownloader(statusFilename="/tmp/pubcrawler.status")
    indexFilename = os.path.join(downloadDir, publisher+".index")

    for metaInfo in readMetaInfos(indexFilename):
        metaInfoDict = metaInfo._asdict()
        pmid = metaInfo.pmid
        textInfo = downloader.downloadFulltext(pmid)

        if textInfo==None:
            logging.debug("No information found, skipping this article")
            continue
        count = 0
        metaInfoDict["fulltextUrl"] = textInfo.baseUrl
        for url, fileType, isSuppdata, httpData in textInfo.getData():
            origFilename = url.split("/")[-1]
            metaInfoDict = fulltextWriter.writeFile(metaInfoDict, str(count), fileType, origFilename, httpData, isSuppdata)
            fulltextWriter.writeMetaInfoDict(metaInfoDict)

def parsePubmedArticle(xmlParser):
    """
    >>> xml = PubmedTestDoc()
    >>> list(parsePubmedArticle(xml))
    [{'pii': 'awq087', 'doi': '10.1093/brain/awq087', 'day-pubmed': '31', 'month-pubmed': '05', 'title': 'Tyrosine hydroxylase deficiency: a treatable disorder of brain catecholamine biosynthesis.', 'vol': '133', 'meshHeadings': 'Age of Onset,Useless Research', 'abstract': 'An infantile onset, progressive, hypokinetic-rigid syndrome with dystonia (type A), and a complex encephalopathy with neonatal onset (type B). Decreased cerebrospinal fluid concentrations of homovanillic acid and c.698G>A and c.707T>C mutations. Carriership of at least one promotor mutation, however, apparently predicts type A tyrosine hydroxylase deficiency. Most patients with tyrosine hydroxylase deficiency can be successfully treated with l-dopa.', 'authors': 'Willemsen MA, Verbeek MM', 'pmcId': '', 'affiliation': 'Radboud University Nijmegen Medical Centre, Donders Institute for Brain, Cognition and Behaviour, Department of Paediatric Neurology (820 IKNC), PO Box 9101, 6500 HB Nijmegen, The Netherlands. m.willemsen@cukz.umcn.nl', 'year-pubmed': '2010', 'articleType': 'research-article', 'pubmedArticleTypes': "Journal Article,Research Support, Non-U.S. Gov't", 'year': '2010', 'pmid': '20430833', 'eIssn': '1460-2156', 'issue': 'Pt 6', 'journal': 'Brain : a journal of neurology', 'fulltextUrl': 'http://www.ncbi.nlm.nih.gov/pmc/articles/PMC/', 'printIssn': '1460-1234'}]
    """
    data = {}
    medlineData           = xmlParser.getXmlFirst("MedlineCitation")
    data["pmid"]          = medlineData.getTextFirst("PMID")
    data["id"]            = "PMID"+data["pmid"]
    data["fulltextUrl"]   = "http://www.ncbi.nlm.nih.gov/pubmed/%s" % data["pmid"]
    logging.log(5, "PMID %s" % data["pmid"])
    data["year-pubmed"]   = medlineData.getTextFirst("DateCreated/Year")
    data["month-pubmed"]  = medlineData.getTextFirst("DateCreated/Month")
    data["day-pubmed"]    = medlineData.getTextFirst("DateCreated/Day")

    artTree               = medlineData.getXmlFirst("Article")
    data["title"]         = artTree.getTextFirst("ArticleTitle")
    data["abstract"]      = artTree.getTextFirst("Abstract/AbstractText")
    data["affiliation"]   = artTree.getTextFirst("Affiliation")
    
    journalTree = artTree.getXmlFirst("Journal")
    data["eIssn"]       = journalTree.getTextFirst("ISSN", reqAttrDict={"IssnType": 'Electronic'})
    data["printIssn"]   = journalTree.getTextFirst("ISSN", reqAttrDict={"IssnType": 'Print'})
    data["vol"]         = journalTree.getTextFirst("JournalIssue/Volume")
    data["issue"]       = journalTree.getTextFirst("JournalIssue/Issue")
    data["year"]        = journalTree.getTextFirst("JournalIssue/PubDate/Year")
    data["journal"]     = journalTree.getTextFirst("Title")

    authorList  = artTree.getXmlFirst("AuthorList")
    lastNames   = []
    initialList = []
    if authorList!=None:
        authorTrees = authorList.getXmlAll("Author")
        for authorTree in authorTrees:
            lastName = authorTree.getTextFirst("LastName", default="")
            if lastName=="":
                lastName = authorTree.getTextFirst("CollectiveName", default="")
            lastNames.append(lastName)

            initials = authorTree.getTextFirst("Initials", default="")
            initialList.append(initials)

    authors = [lastNames[i]+" "+initialList[i] for i in range(0, min(len(lastNames), len(initialList)))]
    data["authors"]=", ".join(authors)

    articleTypeList = artTree.getTextAll("PublicationTypeList/PublicationType")
    articleTypesString  = ",".join(articleTypeList)

    articleType="research-article"

    if "Review" in articleTypeList:
       articleType = "review"
    if "letter" in articleTypeList:
       articleType = "research-article"

    noResearchArticleTags = ["Bibliography", "Biography", 
        "Case Reports", "Webcasts",  
        "Dictionary", "Directory",
        "Editorial", "Festschrift",
        "Patient Education Handout", "Periodical Index", 
        "Portraits", "Published Erratum", "Scientific Integrity Review"
        "Congresses"]

    for noResearchArticleTag in noResearchArticleTags:
        if noResearchArticleTag in articleTypeList:
            articleType = "other"

    data["articleType"]        = articleType
    data["pubmedArticleTypes"] = articleTypesString

    logging.log(5, "pubmedArticleTypes %s, articleType %s" % (articleTypesString, articleType))

    data["doi"]           = xmlParser.getTextFirst("PubmedData/ArticleIdList/ArticleId", reqAttrDict = {"IdType" : 'doi'}, default="")
    data["pmcId"]         = xmlParser.getTextFirst("PubmedData/ArticleIdList/ArticleId", reqAttrDict = {"IdType" : 'pmc'}, default="").replace("PMC", "")
    data["pii"]           = xmlParser.getTextFirst("PubmedData/ArticleIdList/ArticleId", reqAttrDict = {"IdType" : 'pii'}, default="")
    data["fulltextUrl"]   = "http://www.ncbi.nlm.nih.gov/pmc/articles/PMC%s/" % data["pmcId"]

    meshDescriptors = []
    meshHeadingList       = medlineData.getXmlFirst("MeshHeadingList")
    if meshHeadingList:
        for meshHeadingDescriptor in meshHeadingList.getTextAll("MeshHeading/DescriptorName"):
            meshDescriptors.append(meshHeadingDescriptor)

    data["meshHeadings"] = ",".join(meshDescriptors)

    return data

def parsePubmedArticleSet(xml):
    """
    Parse pubmed xml format and yield as dictionary, see parsePubmedArticle
    """
    logging.debug("Parsing pubmed file")
    xp       = maxXml.XmlParser(string=xml)
    for pubmedArticle in xp.getXmlAll("PubmedArticle"):
        yield parsePubmedArticle(pubmedArticle)

def ncbiEFetchGenerator(ids, dbName="pubmed", tool="", email="maximilianh@gmail.com", debug=False):
    """
    retrieve records in xml format from ncbi, and yield as dictionaries
    >> print len(list(ncbiEFetchGenerator(["9322214"], debug=True)))
    1
    """
    idsLeft=list(set(ids))
    retmax = 500

    downloadCount=0
    while len(idsLeft)!=0:
        #logging.debug( "%s" % (str(idsLeft)))
        retStart = downloadCount
        downloadIds = idsLeft[:min(retmax, len(idsLeft))]
        logging.info("Getting data on %d PMIDs from NCBI, %d PMIDs left to download" % (len(downloadIds), len(idsLeft)))
        url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=%s&tool=%s&email=%s&retmode=xml&rettype=medline&id=%s" % (dbName, tool, email, ",".join(downloadIds))
        logging.debug("Getting %s" % url)
        xml = urllib2.urlopen(url).read()
        idsLeft= list(set(idsLeft).difference(downloadIds))
        for pubmedData in parsePubmedArticleSet(xml):
            downloadCount+=1
            yield pubmedData

def ncbiESearch(query, dbName="pubmed", tool="", email="maximilianh@gmail.com", debug=False):
    """ retrieve pmids for query, returns list of ids 
    >> len(set(ncbiESearch("human genome", debug=True)))
    66009
    """
    if debug:
        logging.getLogger().setLevel(5)
        logging.debug("debug mode activated")

    retmax=100000
    addString=""
    query = urllib.quote(query)

    idsLeft = None
    allPmids = []

    while idsLeft>0 or idsLeft==None:
        logging.debug( "%s, %d" % (str(idsLeft), len(allPmids)))
        url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=%s&tool=retrPubmed&tool=%s&email=%s&term=%s&retstart=%d&retmax=%d%s' % (dbName, tool, email, query, len(allPmids), retmax, addString)
        req = urllib2.Request(url)
        #req.add_header('User-Agent', 'Mozilla/5.0 (X11; U; Linux i686; en-US; rv:1.7.8) Gecko/20050524 Fedora/1.5 Firefox/1.5')
        html = urllib2.urlopen(req)
        logging.debug("Getting "+url+"\n")
        
        pmids = []
        for line in html:
            if line.find("<Count>")!=-1:
                if idsLeft==None: # only use the first "count" line
                    fs = re.split("<.{0,1}Count>", line)
                    idsLeft = int(fs[1])
            if line.find("<Id>")!=-1:
                pmid = line.strip().replace("<Id>","").replace("</Id>", "")
                pmids.append(pmid)
        idsLeft -= len(pmids)
        allPmids.extend(pmids)

    return allPmids


def nxmlHasBody(inData):
    """ try to find out if a PMC xml file has some text in it or if 
        it's just scanned pages """
    try:
        root = etreeFromXml(inData)
        bodies = root.xpath("body")
        scans = root.xpath("body/supplementary-material[@content-type='scanned-pages']")
        if len(bodies)>0 and len(scans)==0:
            logging.debug("Found body tag, no scanned pages within it, seems to contain normal fulltext")
            return True
        else:
            logging.debug("No body tag or only scanned pages: No fulltext")
            return False
    except IOError:
        logging.error("IOError while searching for body tag in xml file")
        return False
    
def stripXmlTags(inData, isNxmlFormat=False, isElsevier=False):
    """ read inFile, strip all XML tags, and write to outFile """

    # do not process PMC files without a body tag
    if isNxmlFormat and not nxmlHasBody(inData):
        return None

    try:
        root = etreeFromXml(inData)
        if isElsevier:
            asciiData = treeToAscii_Elsevier(root)
        else:
            asciiData = treeToAsciiText(root)

        return asciiData
    except IOError:
        logging.error("IOError exception while converting xml to text")
        return None
    except SyntaxError:
        logging.error("Syntax Error while converting xml to text")
        return None

class TcpArticleClient:
    def __init__(self, hostname, port):
        self.socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.socket.connect((hostname, int(port)))
        self.readFile=self.socket.makefile('r',0)

    def _sendCleanPair(self, key, value):
        if value=="":
            return
        key   = key.replace("\t", " ").replace("\n", " ")
        value = value.replace("\t", " ").replace("\n", " ")
        line = key+"\t"+value+"\n"
        logging.debug("Sending: %s" % line)
        self.socket.send(line)

    def sendDocument(self, metaInfo, filename, text, isSupp):

        try:
            for key, value in metaInfo._asdict().iteritems():
                value = str(value)
                key   = str(key)
                if key=="metaInfo.abstract":
                    key="text.abstract"
                elif key=="metaInfo.title":
                    key="text.title"
                self._sendCleanPair("metaInfo."+key, value)
                
            self._sendCleanPair("filename", filename)
            self._sendCleanPair("text.body", text)
            self._sendCleanPair("isSupp", str(isSupp))
            self._sendCleanPair("end", "record")
        except KeyboardInterrupt:
            logging.info("Sending quit command to java process")
            self.sendQuit()
            time.sleep(1)

    def receiveResults(self):
        #data = self.socket.recv(1024)
        logging.debug("waiting to receive data")
        lines = []
        while True:
            line = self.readFile.readline()
            line = line.strip()
            logging.debug("received line: "+line)
            if line.lower().startswith("end"):
                logging.debug("end recognized")
                break 
            lines.append(line)
        logging.debug("finished receiving")
        return lines

    #def loadJavaClass(self, className):
        #self._sendCleanPair("loadClass", className)

    def sendQuit(self):
        self.socket.send("quit\tquit")
    
def startJavaProcess(ftHomeDir, className, classParams, debugMode=False):
    """ start a java process based on the jar files in ftHomeDir 
    and return the pid"""
    paramString=""
    if classParams:
        paramString = "--params %s" % classParams

    jarFiles = [os.path.join(ftHomeDir, x) for x in os.listdir(ftHomeDir) if x.endswith(".jar")]
    classPath = ":".join(jarFiles)

    cmdLine = "java -cp %s ac.uk.man.fttools.FtToolsServer --annotator %s %s" % (classPath,  className, paramString)
    if debugMode==True:
        logging.info("Debug level set for java class")
        cmdLine+=" -d"

    logging.info("Starting up java vm with command %s" % cmdLine)
    javaPid = subprocess.Popen(cmdLine, shell=True).pid
    logging.debug("Waiting for 1 second")
    time.sleep(1)
    return javaPid

def remoteResultsGenerator(asciiDir, publisher, nodeCount, nodeId, textType, testId, hostname, portNo, dots):
    """ send documents from local repo to remote server """
    logging.info("Connecting to host %s, port %d" % (hostname, portNo))
    tcpClient = TcpArticleClient(hostname, portNo)

    logging.debug("Sending documents")
    count = 0
    for metaInfo, isSupp, filename, text in generateText(asciiDir, publisher, nodeCount, nodeId, textType, testId=testId):
        tcpClient.sendDocument(metaInfo, filename, text, isSupp)
        lines = tcpClient.receiveResults()
        count+=1
        if count%dots==0:
            logging.info("Processed %d files" % count)
        yield metaInfo, lines

nonLetterRe = re.compile('\W', re.U) # non-alphanumeric character in unicode set

def createUniqueIds(metaInfoFilenameList):
    """ create unique ids for each paper in format 1stauthorYear(letters)
    and return as a list of (numId, uniqueId) """

    if sys.modules.get("unidecode", None)==None:
        logging.error("The unidecode.py module could not be loaded.")
        logging.error("Please install it. Instructions are on http://pypi.python.org/pypi/Unidecode/0.04.1")
        exit(1)

    uniqueIdCounts={}
    uniqueIds = []

    for filename in metaInfoFilenameList:
        for metaInfo in readMetaInfos(filename):
            authors = metaInfo.authors
            firstAuthor = authors.split(" ")[0]
            year = metaInfo.year

            if len(firstAuthor)==0:
                uniqueId = metaInfo.numId
            else:
                firstAuthor = nonLetterRe.sub("", firstAuthor) # remove nonalphanumeric characters
                firstAuthor = unidecode.unidecode(firstAuthor) # substitute accents
                # test if this works with pmcid 1936429
                stem = firstAuthor+year
                stemCount = uniqueIdCounts.get(stem, 0)
                if stemCount == 0:
                    uniqueIdCounts[stem]=1
                    suffix=""
                else:
                    count = uniqueIdCounts[stem]
                    uniqueIdCounts[stem]+=1 
                    suffix = util.baseN(count) 
                uniqueId = stem+suffix
            uniqueIds.append((metaInfo.numId, uniqueId))

    return uniqueIds
    

def createIndex_Genetics(indexFile):
    if os.path.isfile(indexFile):
        os.rename(indexFile, indexFile+".bak")

    pmids = html.ncbiESearch('"Genetics"[journal]', tool=PROGNAME, email=email)
    open(indexFile, "w").write("\n".join(pmids))
    logging.info("Wrote %d PMIDs to file %s" % (len(pmids), indexFile))

def createIndex_Elsevier(elsevierDir, indexFile):
    logging.info("Searching for .xml files in all subdirs of %s" % elsevierDir)
    inFiles = util.findSubdirFiles(elsevierDir, ".xml")
    open(indexFile, "w").write("\n".join(inFiles))

def createIndex_Pmc(indexFilename, localCacheDir):
    localFileList = os.path.join(localCacheDir, "file_list.txt")
    if os.path.isfile(localFileList):
        logging.info("Copying local pmc ftp index %s to %s" % (localFileList, indexFilename))
        shutil.copyfile(localFileList, indexFilename)
    else:
        url = "ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/file_list.txt"
        logging.info("Downloading %s to %s" % (url, indexFilename))
        downloadFile(url, indexFilename)

def strip_namespace_inplace(etree, namespace=None,remove_from_attr=True):
    """ Takes a parsed ET structure and does an in-place removal of all namespaces,
        or removes a specific namespacem (by its URL).
        
        Can make node searches simpler in structures with unpredictable namespaces
        and in content given to be non-mixed.

        By default does so for node names as well as attribute names.       
        (doesn't remove the namespace definitions, but apparently
         ElementTree serialization omits any that are unused)

        Note that for attributes that are unique only because of namespace,
        this may attributes to be overwritten. 
        For example: <e p:at="bar" at="quu">   would become: <e at="bar">

        I don't think I've seen any XML where this matters, though.
    """
    if namespace==None: # all namespaces                               
        for elem in etree.getiterator():
            tagname = elem.tag
            if not isinstance(elem.tag, basestring):
                continue
            if tagname[0]=='{':
                elem.tag = tagname[ tagname.index('}',1)+1:]

            if remove_from_attr:
                to_delete=[]
                to_set={}
                for attr_name in elem.attrib:
                    if attr_name[0]=='{':
                        old_val = elem.attrib[attr_name]
                        to_delete.append(attr_name)
                        attr_name = attr_name[attr_name.index('}',1)+1:]
                        to_set[attr_name] = old_val
                for key in to_delete:
                    elem.attrib.pop(key)
                elem.attrib.update(to_set)

    else: # asked to remove specific namespace.
        ns = '{%s}' % namespace
        nsl = len(ns)
        for elem in etree.getiterator():
            if elem.tag.startswith(ns):
                elem.tag = elem.tag[nsl:]

            if remove_from_attr:
                to_delete=[]
                to_set={}
                for attr_name in elem.attrib:
                    if attr_name.startswith(ns):
                        old_val = elem.attrib[attr_name]
                        to_delete.append(attr_name)
                        attr_name = attr_name[nsl:]
                        to_set[attr_name] = old_val
                for key in to_delete:
                    elem.attrib.pop(key)
                elem.attrib.update(to_set)

def treeToAscii_Elsevier(tree):
    """ try to convert an elsevier XML file to normal ascii text """
    asciiText = ""
    dp = tree.find("document-properties")
    if dp!=None:
        rawTextEl = dp.find("raw-text")
        if rawTextEl!=None:
            asciiText            = rawTextEl.text

    # search for main article tag
    for articleTag in ELSEVIER_ARTICLE_TAGS:
        articleEl = tree.find(articleTag)
        if articleEl is not None:
            break

    if articleEl is None:
        logging.error("Invalid article tag")
        asciiText = ""
    else:
        if asciiText=="":
            asciiText           = treeToAsciiText(articleEl)
    return asciiText

def parseOneElsevierFile(string, data):
    """
    use elementTree to parse Elsevier Consys XML to the metaData dictionary
    """
    def findText(tree, string):
        el = tree.find(string)
        if el!=None:
            return el.text
        else:
            return ""

    hasFulltext = False

    tree   = etreeFromXml(string)

    desc =  tree.find("RDF/Description")
    if desc==None:
        logging.warn("Uppercase RDF/Description not found.")
        desc =  tree.find("rdf/description")

    # PARSE META INFORMATION
    data["source"]          = "elsevier"
    data["title"]           = findText(desc, "title")
    data["articleType"]     = findText(desc, "aggregationtype")
    data["doi"]             = findText(desc, "doi")
    data["id"]              = data["doi"]
    data["journal"]         = findText(desc, "publicationname")
    data["fulltextUrl"]     = findText(desc, "url")
    data["printIssn"]       = findText(desc, "issn")
    data["page"]            = findText(desc, "startingpage")
    data["issue"]           = findText(desc, "number")
    data["year"]            = findText(desc, "coverdisplaydate")
    data["vol"]             = findText(desc, "volume")

    # authors: first try to get from description
    creator = desc.find("creator")
    authors = []
    if creator is not None:
        cSeq = creator.find("seq")
        if cSeq is not None:
            for liEl in cSeq.iterfind("li"):
                if liEl.text!=None:
                    authors.append(liEl.text)
    data["authors"]="; ".join(authors)

    # try to get year and raw text from properties
    dp = tree.find("document-properties")
    dpYear = dp.find("year-first")
    if dpYear is not None:
        data["year"]            = dpYear.text
    if dp.find("raw-text") is not None:
        hasFulltext=True

    # search for main article tag
    for articleTag in ELSEVIER_ARTICLE_TAGS:
        articleEl = tree.find(articleTag)
        if articleEl is not None:
            break

    # PARSE MAIN ARTICLE
    if articleEl is None:
        logging.warn("No article tag")
    else:
        headEl = articleEl.find("head")
        # overwrite author infos
        if headEl is not None:
            authorGroup = headEl.find("author-group")
            if authorGroup is not None:
                authorNames = []
                for authorEl in authorGroup.iterfind("author"):
                    firstName = authorEl.findtext("given-name")
                    if firstName is None:
                        firstName=""
                    famName   = authorEl.findtext("surname")
                    if famName is None:
                        famName = ""
                    authorNames.append(famName+", "+firstName)
                data["authors"]="; ".join(authorNames)
            
            data["title"] = findText(headEl, "title")

            abstractEl = headEl.find("abstract")
            if abstractEl is not None:
                data["abstract"] = etree.tostring(abstractEl, method="text", encoding="utf8")

            # only extract XML if there is a body tag
            bodyEl = articleEl.find("body")
            if bodyEl is not None:
                hasFulltext=True

    cleanMetaDict = {}

    # fix unicode bug 
    for key, val in data.iteritems():
        if type(val) is not types.UnicodeType:
            val = val.decode("latin1")
        val = val.replace("\t", " ")
        cleanMetaDict[key] = val

    return hasFulltext, cleanMetaDict, asciiText 

def download_Genetics(downloadDir, metaInfoFilename, indexFilename, minNumId, testId):
    pmids = set()

    if testId:
        articleParser.downloadOneGeneticsArticle(testId, downloadDir, metaInfoFilename)
    else:
        oldPmids=set()
        internalIds = set()
        # parse oldpmids and internal ids from metainfo
        if os.path.isfile(metaInfoFilename):
            metaInfoFile    = open(metaInfoFilename)
            firstLineFields = metaInfoFile.readline().strip("\n").strip("#").split("\t")
            pmidIndex       = firstLineFields.index("pmid")
            internalIdIndex = firstLineFields.index("numId")
            for line in metaInfoFile:
                parts = line.split("\t")
                oldPmids.add(parts[pmidIndex])
                internalId = int(parts[internalIdIndex])
                internalIds.add(internalId)
            metaInfoFile.close()

        maxInternalId = max(internalIds)
        logging.debug("Maximum numerical internal id is %d"% maxInternalId)
        # read new ids and remove oldpmids
        pmids = [l.strip() for l in open(indexFilename).readlines()]
        logging.debug("id count before old id removal: %d"% len(pmids))
        pmids = set(pmids).difference(oldPmids) # do not download already downloaded stuff
        logging.debug("id count after old id removal: %d"% len(pmids))
        sleepSeconds = maxConfig.getFloat(publisher, "waitSeconds", 5)
        minInternalId = minNumId

        if testId!=None:
            downloadOneGeneticsArticle(testId, outDir, metaInfoFilename, minNumId)
        else:
            for pmid in pmids:
                try:
                    downloadOneGeneticsArticle(pmid, outDir, metaInfoFilename, minNumId)
                except urllib2.HTTPError, ex:
                    if ex.code==404:
                        logging.error("HTTP Error 404 (unknown pmid), skipping article")
                except:
                    logging.error("Error %s occured during download, pausing for ten minutes" % str(sys.exc_info()))
                    time.sleep(float(600))
                minNumId+=1
                logging.debug("Waiting for %s seconds..." % str(sleepSeconds))
                time.sleep(float(sleepSeconds))

def download_Elsevier(downloadDir, indexFilename, metaInfoFilename, nodeCount, nodeId, testId, minNumId, ignoreNoFulltext, elsevierDir): 
    """
    import files from indexFilename into downloadDir and asciiDir at the
    same time. skip files that are already in metaInfoFilename. 
    """

    def importOneElsevier(filename, metaInfoFilename, downloadDir, numId):
        origFilename=util.relpath(filename,elsevierDir)
        metaInfo  = createEmptyMetaDict(source = "elsevier", numId = numId, origFile=origFilename)

        xmlString = codecs.open(filename, encoding = "utf8").read()
        if len(xmlString)==0:
            logging.debug("Skipping file %s: zero file length" % filename)
            return

        hasFulltext, metaInfo, asciiData  = parseOneElsevierFile(xmlString, metaInfo)

        basename = os.path.splitext(os.path.basename(filename))[0]

        if not hasFulltext and ignoreNoFulltext:
            logging.debug("Skipping file %s, no fulltext found and ignoreNoFulltext-flag is set" % filename)
            return
        else:
            localXmlBasename = generate2LevelPath(downloadDir, basename)

            if hasFulltext:
                localName = localXmlBasename + ".xml"
                logging.debug("Writing XML data for article %s to %s" % (filename, localName))
                codecs.open(localName,"w", encoding="utf8").write(xmlString)
                metaInfo["fulltextFileXml"]=util.relpath(localName, downloadDir)
            #if asciiData!=None:
                #localAsciiBasename=util.relpath(localXmlBasename, downloadDir)
                #localAsciiFilename=os.path.join(asciiDir, localAsciiBasename)+".txt"
                #makeAllSubdirs(os.path.dirname(localAsciiFilename))

                #logging.debug("Writing ASCII data for article %s to %s" % (filename, localAsciiFilename))
                #codecs.open(localAsciiFilename,"w", encoding="utf8").write(asciiData)
                #metaInfo["fulltextFile"]=util.relpath(localAsciiFilename,asciiDir )

            appendMetaInfo(metaInfoFilename, metaInfo) 

    def processLine(line, testId, numId, skipIds):
        filename = line.strip("\n")
        relFilename = util.relpath(filename, elsevierDir)
        if testId != None and relFilename!=testId:
            return
        relFilename = util.relpath(filename, elsevierDir)
        if testId == None and relFilename in skipIds:
            logging.debug("Skipping %s, already processed" % (filename))
            return
        logging.debug("Parsing %s" % filename)
        importOneElsevier(filename, metaInfoFilename, downloadDir, numId, asciiDir)

    # FUNCTION START
    skipIds, numId = readOldMetaFile(metaInfoFilename, minNumId, idField="origFile")
    logging.info("Parsing %s to find files to process" % indexFilename)
    indexFile = open(indexFilename)
    runFunctionOnTextfile(indexFile, processLine, minNumId, nodeCount, nodeId, testId, skipIds)

class FileArticleReader:
    """ a class to read files from individual textfiles, referenced by a metaInfo file """
    def __init__(self, asciiDir, dataset):
        metaInfoFilename = os.path.join(asciiDir, dataset+".metaInfo")
        self.metaInfoFilename = metaInfoFilename
        self.asciiDir = asciiDir

    def config(self, textType=None, nodeCount=None, nodeId=None, testId=None):
        self.textType=textType
        self.nodeCount=nodeCount
        self.nodeId=nodeId
        self.testId=testId

    def iterate(self):
        """ returns (metaInfo, text) tuples """
        for metaInfo in readMetaInfos(self.metaInfoFilename, nodeId=self.nodeId, nodeCount=self.nodeCount):
            if self.testId!=None and metaInfo.id!=self.testId:
                continue
            fulltextFiles = []
            if metaInfo.fulltextFile!="":
                fulltextFiles.append(metaInfo.fulltextFile)
            if metaInfo.fulltextFilePdf!="":
                fulltextFiles.append(metaInfo.fulltextFilePdf)
            if metaInfo.fulltextFileXml!="":
                fulltextFiles.append(metaInfo.fulltextFileXml)

            if self.textType != "all":
                if len(fulltextFiles)==0:
                    logging.info("No fulltext file found, skipping article."+metaInfo.id)
                    continue
                else:
                    fulltextFiles = [fulltextFiles.pop(0)]
                    logging.debug("processing only best file: %s" % fulltextFiles[0])
                
            for fulltextFilename in fulltextFiles:
                relFulltextFilename = fulltextFilename
                fulltextFilename = os.path.join(self.asciiDir, fulltextFilename)
                if os.path.isfile(fulltextFilename):
                    fulltextString = open(fulltextFilename).read()
                else:
                    logging.debug("fulltext file referenced by metaInfo-file does not exist: %s, skipping this file" % fulltextFilename)
                    continue
                yield metaInfo, fulltextString
            
            if metaInfo.suppFiles!="" and (self.textType in ["all","bestsupp"]):
                for relSuppFilename in metaInfo.suppFiles.split(","):
                    suppFilename = os.path.join(self.asciiDir, relSuppFilename)
                    if os.path.isfile(suppFilename):
                        suppString = open(suppFilename).read()
                    else:
                        logging.log(5, "supplementary file referenced by metaInfo-file does not exist: %s, skipping this file" % suppFilename)
                        suppString=""
                    yield metaInfo, suppString

class HttpRunner:
    """ runs a Restful webservice at url with given parametres.
    
    parameters is a dictionary with key=value entries to pass via http
    parameters can reference metaInfo fields like $numId 
    parameters must include the text as a field with value $text.
    """

    def __init__(self, url, parameters, headers, fields=None):
        logging.debug("Setting up algorithm at %s with parameters %s" % (url, str(parameters)))
        self.url = url
        self.paramTemplate = parameters
        if fields!=None:
            self.fields    = [int(x) for x in fields.split(",")]
        else:
            self.fields= None

        self.headers            = headers.split(",")
        self.RowType       = namedtuple.namedtuple("algRecord", self.headers)

    def getHeaders(self):
        return self.headers

    def run(self, metaInfo, text):
        httpParams = {}
        for key, value in self.paramTemplate.iteritems():
            if value=="$text":
                httpParams[key] = text
            elif value.startswith("$"):
                httpParams[key] = metaInfo._asdict()[value.strip("$")]
            else:
                httpParams[key] = value

        logging.log(5, "HTTP POST to %s with parameters %s" % (self.url, str(httpParams)))
        dataEncoded = urllib.urlencode (httpParams)
        req = urllib2.Request(self.url, dataEncoded)
        response = urllib2.urlopen(req)
        data =  response.read()

        lines = data.split("\n")
        for line in lines:
            logging.log(5, "Received data: %s" % line)
  
        rows = []
        for line in lines:
            if "<message" in line or "<error" in line:
                logging.debug("Got message from url: " + line)
                continue
            if len(line)==0:
                continue

            columns = line.split("\t")

            if self.fields!=None:
                newRow = []
                for i in self.fields:
                    newRow.append(columns[i])
                columns = newRow
            
            row = self.RowType(*columns)
            yield row

class FulltextFileWriter:
    """ writes metaInfo and binary fulltext articles to file system """
    def __init__(self, baseDir, datasetName):
        self.metaInfoFilename = os.path.join(baseDir, datasetName+".metaInfo")
        self.baseDir = baseDir

    def writeFile(self, metaInfoDict, fileIdentifier, fileType, origFilename, data, isSupp):
        """ write a file and add a reference to it to the metaInfo Object, 
        fileIdentifier is used to identify the file in an article, like S1 or full 
        """
        extId = metaInfoDict["id"]
        filename = "%s.%s.%s" % (extId,fileIdentifier,fileType)
        filePath = generate2LevelPath(self.baseDir, filename)
        relFilePath = util.relpath(filePath, self.baseDir)

        if isSupp:
            metaInfoDict["suppFiles"] = metaInfoDict["suppFiles"]+","+relFilePath
            metaInfoDict["origSuppFiles"] = metaInfoDict["origSuppFiles"]+","+origFilename
        else:
           assert(metaInfoDict["fulltextFilePdf"] in [None, ""])
           metaInfoDict["fulltextFilePdf"] = relFilePath

        logging.debug("Writing data to %s" % filePath)
        fileObj  = open(filePath, "w")
        fileObj.write(data)
        fileObj.close()
        return metaInfoDict

    def writeMetaInfoDict(self, metaInfo):
        """ write a dictionary of field => data to metaInfo file """
        appendMetaInfo(self.metaInfoFilename, metaInfo)

if __name__ == "__main__":
    import doctest
    doctest.testmod()

