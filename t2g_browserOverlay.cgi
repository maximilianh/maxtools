#!/usr/bin/python
import cgitb, urlparse

# Max libraries
import html, pmcArticle, util

# ---- FUNCTIONS -----

def getUcscTracks(db, bedUrlBase):
    """ return dict with name => (db, bedUrl, address) """
    data = {}
    recs = util.sql(db, "select genome, ucscDb, ucscDbdb.defaultPos from genomes, ucscDbdb where genomes.ucscDb=ucscDbdb.name order by orderKey;")
    for genome, ucscDb, defaultPos in recs:
        bedUrl = urlparse.urljoin(bedUrlBase, ucscDb+".bed")
        data[genome] = (ucscDb, bedUrl, defaultPos)
    return data

def getEnsemblTracks(db):
    """ return dict with name => (db, address) """
    data = {}
    recs = util.sql(db, "select genome, defaultPos from genomes where host like '%ensembl.org' order by genome;")
    for genome, defaultPos in recs:
        data[genome] = defaultPos
    return data

# --- PROCESSING / MAIN

cgitb.enable()
db = pmcArticle.openDb()

# headers
h = html.htmlWriter("stdout")
h.startCgi()
h.head("Text2Genome: Genome Browser overlay", stylesheet="default.css")
h.startBody(menu=pmcArticle.config.menu)

# output UCSC table

h.h2("Text2Genome tracks on genome browsers")

h.h3("UCSC Genome Browser (custom tracks)")

ucscData = getUcscTracks(db, pmcArticle.config.bedDir)
h.startTable([300,600], ["Genome","Link"])
for genome, genomeData in ucscData.iteritems():
    ucscDb, bedUrl, defaultPos = genomeData
    ucscUrl = html.ucscCustomTrackLink("Link", ucscDb, defaultPos, bedUrl)
    h.tableRow([genome+" (%s)" % ucscDb, ucscUrl])
h.endTable()

# output Ensembl table
h.h3("Ensembl Genome Browser (DAS server)")

pmcArticle.readGenomesFromMysql(db)
ensemblData = getEnsemblTracks(db)
h.startTable([300,600], ["Genome","Link"])
for genome, defaultPos in ensemblData.iteritems():
    ensemblUrl = pmcArticle.ensemblUrl(genome, defaultPos, pmcArticle.config.dasUrl)
    ensemblLink = '<a href="%s">link</a>' % ensemblUrl
    h.tableRow([genome, ensemblLink])
h.endTable()

h.endHtml()

