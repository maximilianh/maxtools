# library for convenience functionc related to html ouput
import sys, codecs, urllib2, urllib, cgi, doctest, logging, re, os
try:
    from BeautifulSoup import BeautifulStoneSoup
except:
    #sys.stderr.write("warning html.py: BeautifulSoup not installed.\n")
    pass

try:
    from lxml import etree
except:
    pass


# *** return only the urls
def aniGeneUrl(geneId):
    return 'http://crfb.univ-mrs.fr/aniseed/molecule-gene.php?name=%s' % geneId
def entrezUrl(accNo):
    return 'http://www.ncbi.nlm.nih.gov/sites/entrez?db=Gene&amp;term=%s' % accNo
def ensGeneUrl(geneId, orgName="Homo_sapiens"):
    orgName=orgName.replace(" ", "_")
    return 'http://sep2009.archive.ensembl.org/%s/geneview?gene=%s;db=core' % (orgName, geneId)
def genbankUrl(accId):
    return 'http://www.ncbi.nlm.nih.gov/sites/entrez?db=nuccore&amp;term=%s' % accId
def ncbiTaxonUrl(taxId):
    return 'http://www.ncbi.nlm.nih.gov/sites/entrez?db=taxonomy&amp;term=%s' % str(taxId)
def pmcFulltextXmlUrl(pmcId):
    return 'http://www.pubmedcentral.nih.gov/oai/oai.cgi?verb=GetRecord&identifier=oai:pubmedcentral.nih.gov:%s&metadataPrefix=pmc' % str(pmcId)

def ensemblAutoAttachUrl(genome, chrom, start, end, ourDasUrl):
    coords = "%s:%d-%d" % (chrom, start, end)
    # example: "http://www.ensembl.org/Homo_sapiens/Location/View?g=ENSG00000012048;contigviewbottom=das:http://www.ensembl.org/das/Homo_sapiens.NCBI36.transcript=labels"
    url = "http://sep2009.archive.ensembl.org/"+genome+"/Location/View?r=%s;contigviewbottom=das:%s=labels" %(coords, ourDasUrl)
    return url
def ucscTrackUrl(hgsid, chrom, start, end, server="http://genome.ucsc.edu"):
    return server+"/cgi-bin/hgTracks?hgsid=%s&position=%s:%s-%s" % (hgsid, chrom, start, end)

def pubmedUrl(pmid):
    return 'http://www.ncbi.nlm.nih.gov/pubmed/%s' % (pmid)
def pmcUrl(pmcId):
    return 'http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=%s' % (pmcId)

# *** return complete links = with <a> tags around and a default text

def pubmedLink(pmid, text=None):
    if not text:
        text=pmid
    return '<a href="http://www.ncbi.nlm.nih.gov/pubmed/%s">%s</a>' % (pmid, text)

def pmcLink(articleId, text=None):
    if not text:
        text=articleId
    artStr = str(articleId)
    return '<a href="http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=%s">%s</a>' % (artStr, text)
def ensGeneLink(geneId, orgName="Homo_sapiens"):
    return '<a href="http://sep2009.archive.ensembl.org/%s/geneview?gene=%s;db=core">%s</a>' % (orgName, geneId, geneId)
def geneCardsLink(symbol):
    return '<a href="http://www.genecards.org/cgi-bin/cardsearch.pl?search=disease&symbols=%s#MINICARDS">%s</a>' % (symbol, symbol)
def aniGeneLink(geneId):
    return '<a href="http://crfb.univ-mrs.fr/aniseed/molecule-gene.php?name=%s">%s</a>' % (geneId, geneId)
def aniInsituPageLink(insituId):
    return '<a href="http://crfb.univ-mrs.fr/aniseed/insitu.php?id=%s">%s</a>' % (insituId, insituId)
def aniISHLink(geneId):
    return '<a href="http://crfb.univ-mrs.fr/aniseed/insitu-result.php?target=%s&BOOLmut=3&BOOLmanip=2&MOLtype=2&Order=DEV_STAGE_ID">%s</a>' % (geneId, geneId)
def genbankLink(accId):
    return '<a href="http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Search&amp;db=Nucleotide&term=%s&amp;doptcmdl=GenBank">%s</a>' %(accId, accId)

def ensGenomeLink(orgName, chrom, start, end):
    coords = '%s:%d-%d' % (chrom, start, end)
    return '<a href="http://www.ensembl.org/%s/Location/View?r=%s">%s</a>' % (orgName, coords, coords)

def ucscGenomeLink(baseUrl, db, pos, hgsid=None, desc=None):
    hgsidStr=""
    if desc==None:
       desc = pos
    if hgsid!=None:
        hgsidStr = "&hgsid=%s" % hgsid
    return '<a href="%s/cgi-bin/hgTracks?db=%s&position=%s%s">%s</a>' % (baseUrl, db, pos, hgsidStr, desc)

def ucscCustomTrackUrl(db, pos, customTrackUrl, server="http://genome.ucsc.edu"):
    return server+"/cgi-bin/hgTracks?db=%s&position=%s&hgt.customText=%s" % (db, pos, customTrackUrl)

def ucscCustomTrackLink(description, db, pos, customTrackUrl, server="http://genome.ucsc.edu"):
    return '<a href="%s">%s</a>' % (ucscCustomTrackUrl(db, pos, customTrackUrl, server), description)

def ucscMafLink(text, baseUrl, db, track, chrom, start, end):
    return '<a href="%s/cgi-bin/hgc?o=%d&t=%d&g=%s&c=%s&l=%d&r=%d&db=%s">%s</a>' % (baseUrl, start, end, track, chrom, start, end, db, text)
def pmcFulltextXmlLink(pmcId):
    return '<a href="%s">PMC%s</a>' % (pmcFulltextXmlUrl(pmcId), str(pmcId))

def ucscHgsid(db, server="http://genome.ucsc.edu"):
    """ get a new hgsid from ucsc """
    """ db is a ucsc database as a string, e.g. hg18"""
    #print("Requesting a new hgsid from UCSC")
    data = urllib2.urlopen(server+"/cgi-bin/hgCustom?db=%s" % db)
    for line in data:
        if line.startswith('<INPUT TYPE=HIDDEN NAME="hgsid" VALUE="') or line.startswith("<INPUT TYPE=HIDDEN NAME='hgsid' VALUE=") :
            line = line.strip()
            line = line.replace('<INPUT TYPE=HIDDEN NAME="hgsid" VALUE="', '')
            line = line.replace('"><TABLE BORDER=0>', '')
            line = line.replace("<INPUT TYPE=HIDDEN NAME='hgsid' VALUE='", '')
            line = line.replace("'>", '')
            #print("New hgsid %s" % line)
            return line
    #sys.stderr.write("error in UCSC web parser, write to maximilianh@gmail.com to get this fixed")
    print("error in UCSC web parser, write to maximilianh@gmail.com to get this fixed")

def ucscUpload(db, data, hgsid=None, server="http://genome.ucsc.edu", name="User track", description="User track", visibility="1", clade="deuterostome", organism="C. intestinalis"):
    """ adds data as a user track, creates and returns session id if hgsid is None, returns None on error """
    """ db is a string like hg18, data is your custom track as one big string (including the newlines)"""
    #log("Uploading %d lines to UCSC" % (data.count("\n")+1))

    if not data.startswith("track"):
        data = 'track name="%s" description="%s" visibility=%s\n' % (name, description, visibility) + data

    if hgsid==None:
        hgsid = ucscHgsid(db, server)

    vars = {}
    vars["hgsid"]=hgsid
    vars["clade"]=clade
    vars["org"]=organism
    vars["db"]=db
    vars["hgct_customText"]=data
    vars["Submit"]="Submit"
    try:
        fullUrl = server+"/cgi-bin/hgCustom/?db="+db
        html = urllib2.urlopen(fullUrl, urllib.urlencode(vars))
    except urllib2.HTTPError:
        print("WARNING: Could not upload track into UCSC genome browser, links to UCSC will therefore not display any matches")
        return None
    html = html.readlines()
    for l in html:
        if l.find("Manage Custom Tracks")!=-1:
            return hgsid
    print "ERROR: Could not upload custom track into UCSC server at %s<br>\n" % server
    print "Offending data was:<br>\n"
    print data
    print "Error Message was:<br>\n"
    print "\n".join(html)
    return None

class htmlWriter:
    def __init__(self,fname):
        if fname=="stdout":
            #codecs.getwriter('utf8')(sys.stdout.buffer)
            self.f = sys.stdout
        else:
            self.f = open(fname, "w")

    def startCgi(self, contentType="text/html; charset=utf-8", addLines=None):
        self.writeLn ("Content-type: %s" % contentType)
        if addLines:
            for l in addLines:
                self.writeLn(l)
        self.writeLn("")

    def write(self, text):
        #self.f.write(text.encode("latin1", 'replace'))
        self.f.write(text.encode("utf", 'replace'))

    def writeLn(self, str):
        if str==None:
            str="None"
        self.write(str)
        self.write("\n")

    def head(self, title, stylesheet=None, analyticsAccount=None):
        #<meta name="keywords" content="Casey Bergman, Casey M. Bergman, Drosophila, noncoding, cis-regulatory, yeast, genome annotation, text mining, transposable element, molecular evolution, TE, genome evolution" />
        #<meta name="description" content="The lab page of Casey M. Bergman" />

        self.f.write("""
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
<html>
<head>
<meta http-equiv="content-type" content="text/html; charset=utf-8" >
<title>"""+title+"""</title>""")
        if stylesheet!=None:
            self.f.write("""<link href="%s" rel="stylesheet" type="text/css" >\n""" % stylesheet)

        if analyticsAccount!=None:
            self.f.write(""" 
           <script type="text/javascript">

              var _gaq = _gaq || [];
              _gaq.push(['_setAccount', '%s']);
              _gaq.push(['_trackPageview']);

              (function() {
                var ga = document.createElement('script'); ga.type = 'text/javascript'; ga.async = true;
                ga.src = ('https:' == document.location.protocol ? 'https://ssl' : 'http://www') + '.google-analytics.com/ga.js';
                var s = document.getElementsByTagName('script')[0]; s.parentNode.insertBefore(ga, s);
              })();

            </script> 
            """ % analyticsAccount)

        self.f.write ("""</head>\n""")


    def endHtml(self):
        self.f.write("\n</body>\n</html>")
        self.f.close()

    def insertMenu(self, menu):
        """ menu: dict with title => list of (url, description)"""
        self.writeLn("")
        self.writeLn('<div class="urbangreymenu">')
        for title, titleData in menu:
            self.writeLn('<h3 class="headerbar">%s</h3>' % title)
            self.writeLn("<ul>")
            for desc, url in titleData:
                self.writeLn('<li><a href="%s">%s</a></li>' % (url, desc))
            self.writeLn('</ul>')
        self.writeLn("</div>")
        self.writeLn("")

    def startBody(self, title=None, menu=None):
        self.writeLn("<body>")
        if menu:
            self.insertMenu(menu)
        if title:
            self.f.write("""<h2>"""+title+"""</h2>""")

    def img(self, url, alt, width=None, height=None):
        if width:
            attr = 'width="%s"' % str(width)
        if height:
            attr += ' height="%s"' % str(height)
        self.f.write('<img src="%s" alt="%s" %s/>\n' % (url, alt, attr))

    def h2(self, text):
        unicode=('<h2>'+text+'</h2>\n')
        self.write(unicode)

    def h3(self, text):
        unicode=('<h3>'+text+'</h3>\n')
        self.write(unicode)

    def h4(self, text):
        unicode=('<h4>'+text+'</h4>\n')
        self.write(unicode)

    def linkList(self, listlst):
        self.f.write('Index:<ul>\n')
        for pair in listlist:
            desc, name = pair
            self.f.write('    <li>\n')
            self.f.write('    <a href="#%s">%s</a>\n' % (name, desc))
        self.f.write('</ul>\n')

    ### TABLES

    def startTable(self, widths, headers, bgcolor=None, cellspacing=2, cellpadding=3):
        if bgcolor==None:
            self.f.write("""\n<table  border="0" cellspacing="%d"  cellpadding="%d"  >\n""" % (cellspacing, cellpadding))
        else:
            self.f.write('\n<table  border="0" cellspacing="%d"  cellpadding="%d"  bgcolor="%s">\n' % (cellspacing, cellpadding, bgcolor) )
        #rules="cols"frame="hsides"
        if len(widths)>0:
            self.f.write("<colgroup>\n")
            for width in widths:
                self.f.write("  <col width='%s'>\n" % (str(width)))
            self.f.write("</colgroup>\n\n")
        if len(headers)>0:
            self.f.write("  <tr>\n")
            for header in headers:
                self.f.write('    <th align="left">%s</th>\n' % (header))
            self.f.write("  </tr>\n")

    def endTable(self):
        self.f.write("</table>\n")

    def startTd(self):
            self.f.write('<td>\n')

    def endTd(self):
            self.f.write('</td>\n')

    def tableData(self, str, colour=None):
        if colour==None:
            self.f.write('<td>%s</td>\n' % str)
        else:
            self.f.write('<td bgcolor="%s">%s</td>\n' % (colour, str))

    def startTr(self, bgcolor=None):
        if bgcolor==None:
            self.f.write('<tr valign="top">\n')
        else:
            self.f.write('<tr bgcolor="%s" valign="top">\n' % bgcolor)

    def endTr(self):
        self.f.write("</tr>\n")

    def startTd(self):
        self.f.write("<td>\n")

    def endTd(self):
        self.f.write("</td>\n")

    def tableRow(self, cellList, bgColour=None, colorList=None):
        self.startTr(bgColour)

        for i in range(0, len(cellList)):
            if colorList!=None:
                col=colorList[i]
            else:
                col=None
            cell = cellList[i]
            self.tableData(cell, col)
        self.endTr()

    ## TEXT FORMATTING 
    def link(self, url, text):
        self.write('<a href="%s">%s</a>' % (url, text))

    def linkStr(self, url, text):
        return '<a href="%s">%s</a>' % (url, text)

    def anchor(self, name):
        self.f.write('<a name="%s"></a>' % name)

    def h3(self, str, anchor=None):
        if anchor!=None:
            self.f.write('<a name="%s"></a>\n' % anchor)
        self.f.write('<hr>\n<h3>%s</h3>\n' % str)

    def small(self,str):
        self.f.write('<small>%s</small>' % str)

    def b(self, str):
        self.f.write('<b>%s</b>' % str)

    def br(self):
        self.f.write("<br>\n")

    def pre(self, text):
        self.f.write("<pre>%s</pre>\n" % text)

    def p(self):
        self.f.write("<p>\n")

    def centerStart(self):
        self.f.write("<center>\n")

    def centerEnd(self):
        self.f.write("</center>\n")

    def startUl(self):
        self.f.write("<ul>\n")

    def li(self, text):
        self.f.write("<li>%s</li>\n" % text)

    def endUl(self):
        self.f.write("</ul>\n")


    ## FORMS 

    def formStart(self, action, method="get"):
        self.writeLn('<form name="input>" action="%s" method="%s">' % (action, method))

    def formInput(self, type, name, size=None, value=None):
        addStr=""
        if size:
            addStr+='size="%d"' % size
        if value:
            addStr+='value="%s"' % value

        self.writeLn('<input type="%s" name="%s" %s />' % (type, name, addStr))

    def formInputText(self, name, size=None):
        self.formInput("text", name, size)

    def formInputSubmit(self, name):
        self.formInput("submit", name, value=name)

    def formInputReset(self, name):
        self.formInput("reset", name)

    def formEnd(self):
        self.writeLn('</form>')

    ## LINKS TO EXTERNAL RESOURCES
    def ensGeneLink(self, orgName, geneId):
	return ensGeneLink(orgName, geneId)

    def ucscGenomeUrl(self, baseUrl, db, pos):
        return '%s/cgi-bin/hgTracks?db=%s&position=%s' % (baseUrl, db, pos)

    def zfinGeneLink(self, geneId, title=None):
	if title==None:
	    title=""
	else:
	    title = ' title="%s" ' % title
	return '<a href="http://zfin.org/cgi-bin/webdriver?MIval=aa-markerview.apg&OID=%s" %s>%s</a>' % (geneId, title, geneId)

    def geneCardsLink(self, symbol):
	return geneCardsLink(symbol)

    def brainAtlasLink(self, mouseEntrezIds):
	abaItems = ["entrezgeneid:=" + i for i in mouseEntrezIds]
	query = " OR ".join(abaItems)
	ids = ", ".join(mouseEntrezIds)
	return '<a href="http://www.brain-map.org/search.do?findButton.x=1&queryText=%s">%s</a>' % (query, ids)

    def zfinInsituLink(self, geneId, desc=None, title=None):
	if desc==None:
	    desc=geneId
	if title!=None:
	    titleStr=' title="%s" ' % title
	else:
	    titleStr=""
	return '<a href="http://zfin.org/cgi-bin/webdriver?MIval=aa-xpatselect.apg&query_results=true&xpatsel_geneZdbId=%s" %s">%s</a>' % (geneId, titleStr, desc)

    def ghostGeneUrl(self, geneId):
	return 'http://ghost.zool.kyoto-u.ac.jp/cgi-bin3/txtgetr2.cgi?%s' % geneId

    def ghostInSituUrl(self, geneId):
	return 'http://ghost.zool.kyoto-u.ac.jp/cgi-bin3/photoget2.cgi?%s' % geneId

    def aniseedGeneUrl(self, geneId):
	return aniseedGeneUrl(geneId)

    def aniseedInSituUrl(self, geneId):
	return 'http://crfb.univ-mrs.fr/aniseed/insitu-result.php?target=%s&BOOLmut=3&BOOLmanip=2&MOLtype=2&Order=DEV_STAGE_ID' % geneId

    def ensemblMCUrl(self, baseName, compName, basePair, compPair):
        """ return string with html-link to ensembl multicontigview from two genes given genes and organism names """
        urlMask = "http://www.ensembl.org/%s/multicontigview?s1=%s;w=%d;c=%s:%d:1;w1=%d;c1=%s:%d:1;action=%s;id=1"  
        baseSize = (basePair.right.end - basePair.left.start) * 2
        baseChrom = basePair.left.chrom.replace("chr","")
        basePos = basePair.left.start
        compChrom = compPair.left.chrom.replace("chr","")
        compPos = compPair.left.start
        if compPair.left.start - compPair.right.start < 0:
            action="out"
            compSize = (compPair.right.end - compPair.left.start) * 2
        else:
            action="flip"
            compSize = (compPair.left.end - compPair.right.start) * 2
        url = urlMask % (baseName, compName, baseSize, baseChrom, basePos, compSize, compChrom, compPos, action)
        text = "(Multicontigv.)"
        urlStr = '<a href="%s">%s</a>\n' % (url, text)
        return urlStr

## CGI
    def gotCgiVariables(self):
        self.cgiVarsRaw = cgi.FieldStorage()

        self.cgiVars={}
        for var in self.cgiVarsRaw:
            self.cgiVars[var]=self.cgiVarsRaw[var].value

        return len(self.cgiVars)!=0

def HTMLEntitiesToUnicode(text):
    """Converts HTML entities to unicode.  For example '&amp;' becomes '&'."""
    text = unicode(BeautifulStoneSoup(text, convertEntities=BeautifulStoneSoup.ALL_ENTITIES))
    return text

def unicodeToHTMLEntities(text):
    """Converts unicode to HTML entities.  For example '&' becomes '&amp;'."""
    text = cgi.escape(text).encode('ascii', 'xmlcharrefreplace')
    return text

if __name__ == "__main__":
    import doctest
    doctest.testmod()

