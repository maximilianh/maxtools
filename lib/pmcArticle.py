# std python libs
import operator, random, logging

# import t2g config and add lib directory to python module search path
import t2gConfig
# external library, installed on agent (manchester), not on vital-it cluster (lausanne)
try:
    import MySQLdb
except:
    pass # silently ignore error message
    
import tabfile

import util, pmcArticle

# TODO:
# - filtering could be improved: only remove sequences of which the BEST match is in univec ? currently, any match with univec removes the whole sequence

from sys import *

# === FUNCTIONS ==================

def getSeqDataFromMysql(db, pmcId):
    #db = MySQLdb.connect(host, user, pwd, dbase, port)
    pmcObj = PmcArticle()
    res = pmcObj.readFromMysql(db, pmcId)
    if res==None:
        return None
    else:
        return pmcObj

def openDb(**kwargs):
    #db = MySQLdb.connect(t2gConfig.mysqlHost, t2gConfig.mysqlUser, t2gConfig.mysqlPassword, t2gConfig.mysqlDb, t2gConfig.mysqlPort, **kwargs)
    db = MySQLdb.connect(t2gConfig.mysqlHost, t2gConfig.mysqlUser, t2gConfig.mysqlPassword, t2gConfig.mysqlDb, t2gConfig.mysqlPort, read_default_file="/etc/my.cnf", charset="utf8", **kwargs)
    db.set_character_set("utf8")
    return db

def getNumIds(db, geneId, asStrings=False):
    recs = util.sql(db, "select numId from bestGenes where geneId='%s';" % geneId)
    recs = [l[0] for l in recs]
    if asStrings:
        recs = [str(x) for x in recs]
    return recs

def getHumanOrthologs(db, geneId):
    recs = util.sql(db, "select hsGeneId from orthologs where geneId='%s' " % (geneId))
    recs = [l[0] for l in recs]
    return recs

def getIdForSymbol(db, symbol):
    recs = util.sql(db, "select geneId from genes where symbol = %s;", symbol)
    if len(recs)==0:
        return symbol
    else:
        return recs[0][0]

def getOrthologs(db, geneIds):
    geneIds = ["'"+g+"'" for g in geneIds]
    query = "select geneId from orthologs where hsGeneId in (%s);" % ",".join(geneIds)
    recs = util.sql(db, query)
    recs = [l[0] for l in recs]
    return recs

def getNumIdsOrthologs(db, orthologs):
    orthologIds = ["'"+g+"'" for g in orthologs]
    query = "select numId from bestGenes, orthologs where orthologs.geneId in (%s) and bestGenes.geneId=orthologs.geneId;" % ",".join(orthologIds)
    recs = util.sql(db, query)
    recs = [str(l[0]) for l in recs]
    return recs

#def sql(db, sql):
    #if db==None:
        #db = openDb()
#
    #cursor = db.cursor()
    #cursor.execute(sql)
    #data = cursor.fetchall()
    #cursor.close ()
    #db.commit()
    #return data

def bestScoreElements(list, scoreGetter=operator.itemgetter(0), returnCompleteList=False):
    """ in: list of tuples (score, element) will sort by score and return only those elements with a score as good as the best element """
    """ alternative: list of tuples, scoreGetter can retrieve score from tuples, will return tuples where score >= bestScore """
    list.sort(key=scoreGetter, reverse=True)
    bestScore = scoreGetter(list[0])
    if returnCompleteList:
        bestElements = [e for e in list if scoreGetter(e) >= bestScore]
        return bestElements
    else:
        bestElements = [e for x,e in list if x >= bestScore]
        return set(bestElements)


# for speed reasons: global variable in module
# so we are reading this table only once per http request
idToGenome={}
genomeToId = {}
genomeToAssembly = {}

#def readGenomesFromMysql(db):
    #genomesRecs = util.sql(db, "select genome, assembly, genomeId from genomes;")
    #for genome, assembly, id in genomesRecs:
        #id = int(id)
        #idToGenome[int(id)] = genome
        #genomeToId[genome] = id
        #genomeToAssembly[genome] = assembly

def ensemblUrl(genome, genomeId, pos, dasUrl=t2gConfig.dasUrl, baseServer="www.ensembl.org"):
    genome = genome.replace(" ","_")
    ourDasUrl = dasUrl+"/"+t2gConfig.dasUri(genomeId)
    # example: "http://www.ensembl.org/Homo_sapiens/Location/View?g=ENSG00000012048;contigviewbottom=das:http://www.ensembl.org/das/Homo_sapiens.NCBI36.transcript=labels"
    url = "http://%s/"%baseServer+genome+"/Location/View?r=%s;contigviewbottom=das:%s=labels" %(pos, ourDasUrl)
    return url

# -----

class BestTargetPredictor:
    """ a class with methods that predict the best target (=genome or gene) given a (pmid, count)-scoring file (=papers per target, attr of the object) 
    and give the hits (parameters to individual methods). The different methods implement different typs of predictions. """

    # attributes: genomePmidCount = {}
    # problemPmids = 0

    def __init__(self, fname="none"):
        self.readPmidCountFile(fname)

    def readPmidCountFile(self, fname):
        """ load a list of (taxonId, count)-tuples to populate the genomePmidCount """
        self.genomePmidCount = {}
        self.problemPmids = 0
        if fname!="none":
            for l in open(fname):
                if l.startswith("#"):
                    continue
                fs = l.strip().split("\t")
                genomeid, pmidCount = fs
                pmidCount = int(pmidCount)
                genomeid = int(genomeid)
                assert(genomeid not in self.genomePmidCount) # make sure that there are no duplicates in genome list
                self.genomePmidCount[genomeid]=pmidCount

    def readBestGenomeFile(self, fname):
        """ load a list of (pmcid, bestGenomeId) tuples for later use by bestGenome_fromFile """
        self.bestGenome = tabfile.slurpdict(fname)

    def bestGenome_fromFile(self, pmcId, hits):
        return self.bestGenome[pmcId] # if you get an error here, there is a serious problem with your bestGenome file

    def bestGenome_mostSequences(self, pmcId, hits):
        """ search for the genome with the highest number of matching sequences
        using only the best hits per sequence
        if several results are found use genome with best ncbi gene score """
        # iterate over seqid, increase genomeScore for each genome 
        genomeScores = {}
        for seqId, seqHits in hits.iteritems():
            genomes = set([h[0] for h in seqHits])
            for g in genomes:
                genomeScores.setdefault(g, 0)
                genomeScores[g]+=1

        genomeScores = genomeScores.items() # convert to list (hashes are not sortable)
        genomeScores = [(y,x) for x,y in genomeScores] # convert to (score, genome) format

        return genomeScores
        # find genome with highest score and return best one, weighted with NCBI Gene
        #bestScoreGenomes = bestScoreElements(genomeScores)
        #bestGenome, bestWeightGenomes = self.mostProbableGenome(bestGenomes)
        #return bestGenome, bestScoreGenomes, bestWeightGenomes

    def bestGenome_bestHit(self, pmcId, hits):
        """ search for genome with single best match, keep only hits with this score 
        if several hits are found use genome with best ncbi gene score """

        def asFlatList(hits):
            """ return dict seqId => list of tuples as flat list of tuples seqId, genome, chrom, start, end, score """
            list = []
            for seqId, seqHits in hits.iteritems():
                for h in seqHits:
                    genome, chrom, start, end, score, percId = h
                    list.append((seqId, genome, chrom, start,end, score, percId))
            return list

        genomeScores = []
        hits = asFlatList(hits)
        for h in hits:
            genome = h[1]
            score = h[-1]
            genomeScores.append((score, genome))

        # sort list by score and keep all elements with score >= max score of list
        #bestGenomes = bestScoreElements(genomeScores)
        # return most probable genome according to NCBI genes
        #return self.mostProbableGenome(bestGenomes), genomeScores
        return genomeScores

    def bestGenome_bestTotalSum(self, pmcId, hits):
        """ search for highest sum of scores per genome, using the highest score per sequence, returning most probable genome """
   
        genomeScores = []
        for genome, genomeHits in hits.byGenome().iteritems():
            genomeScore = 0
            for (seqId, genome, chrom, start, end, score, percId) in genomeHits:
                genomeScore += score
            genomeScores.append( (genomeScore, genome) )
        #bestGenomes = bestScoreElements(genomeScores)
        #return self.mostProbableGenome(bestGenomes), genomeScores
        return genomeScores

    def weightGenomes(self, genomes):
        """ given a list of genomes, return the one that occurs the most often in NCBI gene , """ 
        assert(len(self.genomePmidCount)!=0) # you have to call readPmidCountFile before calling this function
        genomeScores = []
        for g in genomes:
            score = self.genomePmidCount.get(g,0) # if not in file, use score 0
            genomeScores.append( (score, g) )
        bestElements = bestScoreElements (genomeScores)
        #infoMsg=None

        # cannot decide -> use random element
        if len(bestElements)>1:
            bes = [str(e) for e in bestElements]
            self.problemPmids+=1
            bestElements = list(bestElements)
            random.shuffle(bestElements)
            bestElement = bestElements[0]
        # otherwise take any element
        else:
            bestElement = list(bestElements)[0]

        return (bestElement, bestElements, genomeScores)

# ---------

#class BestGenePredictor:
    #""" get best gene given a list of blast transcript hits and a mysql or textfile with bestGenome information """
#
    #def init(self, fname=None, db=None):
        #""" load text file with gene<tab>weight data, if 2 genes have the same support, the one with the better weight will be chosen"""
        #if fname=="mysql" or fname==None:
            #if db==None:
                #self.db = openDb()
            #else:
                #self.db=db
        #elif "geneWeights" not in self.__dict__ or self.geneWeights==None or len(self.geneWeights)==0:
            #self.geneWeights = tabfile.slurpdict(fname)
#
    #def getGeneWeight(self, gene):
        #""" get gene weight from self.geneWeights or alternatively from mysql table named geneWeights, depending 
            #on how loadGeneWeights was called before"""
        #if "db" in self.__dict__:
            #sqlStr = "select weight from geneWeights where geneId='%s';" % gene
            #weightRecs = sql(self.db, sqlStr)
            #assert(len(weightRecs)<=1)
            #if len(weightRecs)>0:
                #return int(weightRecs[0][0])
            #else:
                #return 0
        #else:
            #return self.geneWeights.get(gene, 0)

def bestGenes(hits, geneWeights):
    """ return the set of predicted gene models from one genome that have >=2 filtered blast matches"""
    # create dict with gene -> sequences that match it
    geneScores = {}
    for tuple in hits:
            seqId, genomeId, gene, start, end, score = tuple
            if "|" in gene:
                fs = gene.split("|")
                if len(fs)==2:
                    gene = fs[1]
            geneScores.setdefault(gene, set())
            geneScores[gene].add(seqId)

    # keep only genes with >= 2 matches
    seqIds = {}
    toDel = set()
    for gene, seqSet in geneScores.iteritems():
        seqSetStr = [str(seq) for seq in seqSet]
        if len(seqSet)<2:
            toDel.add(gene)
        else: 
            # save sequence ids that match 
            seqIds[gene]=",".join(seqSetStr)
    for gene in toDel:
        del geneScores[gene]

    # now inverse strategy: create dict with supporting sequences -> genes
    seqSetToGenes = {}
    for gene, seqSet in geneScores.iteritems():
        seqSetToGenes.setdefault(str(seqSet), set()).add(gene)

    # for each seqSet: get gene weight, sort by weight and take gene with best weight
    #predGenes = SetWithAttributes()
    predGenes = set()
    undecGenes = []# number of seqSeqs that match two genes with exactly the same score

    # if identical best weights, return all of them
    for seqSet, geneSet in seqSetToGenes.iteritems():
        geneScores = []
        #print seqSet
        for gene in geneSet:
            #geneWeight = self.getGeneWeight(gene)
            geneWeight = geneWeights.get(gene, 0)
            #print gene, geneWeight
            geneScores.append( (geneWeight, gene) )
        bestGenes = bestScoreElements(geneScores) 
        #print "best", bestGenes
        if len(bestGenes)>1:
            undecGenes.append(bestGenes)

        predGenes.update(bestGenes)

        supportList = []
        for bestGene in bestGenes:
            for score, gene in geneScores:
                desc = "%s=%d" % (gene, score)
                supportList.append(desc)

    #predGenes.seqIds = scoreDetails
    #predGenes.undecGenes = undecGenes
    return predGenes, undecGenes, seqIds

class IndexedBlastHits(dict):
    """ a hash of blast hits from sequences of one pmcId: sourceSequenceId ==> list of  (genome, chrom, start, end, score) """

    def __init__(self, tuples=[]):
        """ construct an object from a list of tuples from e.g. mysql or textfile """
        """ tuples is a list of tuples (pmcId, genome, seqId, chrom, start, end, score) """
        """ will flag all entries with matching genome == 0 (univec) """
        lastPmcId=None
        for t in tuples:
            pmcId, genome, seqId, chrom, start, end, score= t[:7]
            if len(t)>6:
                percId = t[6]
            else:
                percId = 0
            assert(pmcId==lastPmcId or lastPmcId==None) # a set of tuples MUST HAVE identical pmcId fields
            lastPmcId==pmcId
            genome = int(genome)
            seqId = int(seqId)
            #if genome in t2gConfig.filterGenomeIds:
            self.setdefault(seqId, []).append([genome, chrom, int(start),int(end),int(float(score)), int(percId)])
        # this is just to speed up filter() see below, so we don't have to iterate twice over all features

    def onlyBestHits(self):
        # keep only best hits per sequence and put into new hits obj
        filteredHits = IndexedBlastHits() 
        for seqId in self:
            filteredHits[seqId] = bestScoreElements(self[seqId], scoreGetter=operator.itemgetter(4), returnCompleteList=True)
        return filteredHits

    def filter(self):
        """ keep only BEST hits (highest scores) and non-univec sequences for each target, returns a new IndexedBlastHits-object"""
        # input: dict SeqId -> list of [genome, chrom, start, end, score]
        # remove sequences that match a filter-genome
        #seqIds = set(self.keys()).difference(self.filterSeqs)
        filteredHits = self.onlyBestHits()

        # search and remove sequences that contain best matches for univec 
        filterIds = set()
        for seqId, hits in self.iteritems():
            for h in hits:
                if h[0]==0: # univec taxonid is set to 0 in blastConvert
                    filterIds.add(seqId)

        for id in filterIds:
            del filteredHits[id]
        return filteredHits # return new object

    def toString(self, taxonToName):
        data = []
        for key, values in self.iteritems():
            data.append("   seq No: "+str(key))
            for v in values:
                genome, chrom, start, end, score = v
                line = [taxonToName.get(genome, str(genome)), str(chrom), str(start), str(end), "score="+str(score)]
                data.append("      "+", ".join(line))
        #data.append("SeqNo's to remove: %s" % str(self.filterSeqs))
        return "\n".join(data)

    def chain(self, genome, maxDist):
        """ fuse features into groups that are on same genome/transcript and closer than maxDist bp. 
        Return in a format that resembles 12-field-bed:
        chrom, start, end, pmcId, SequenceIDS(extended name), +, start, end, 0, blockSizeLen, blockSizes, blockStart
        so the name field is split into fields 3 and 4 (zero based)
        return format is list of [chrom, start, end, seqid, score, seqIds, groupId] )
        For genes, we do not create a group id but save a copy of the transcript name instead
        """
        fusedSeqs = set()
        maxInt = 2**32

        startIdx, endIdx, seqIdx, scoreIdx = 0, 1, 2, 3

        groupedFeatures = []
        chromFts = self.byChrom(genome) # use only features from bestGenome
        for chrom, features in chromFts.iteritems():
            features.sort(key=lambda f: f[0]) # sort fts by start pos
            clusters = []
            lastEnd = None

            for f in features:
                start, end, seqid, score = f
                if lastEnd==None or (start - lastEnd) > maxDist: 
                    clusters.append([])
                clusters[-1].append(f)
                lastEnd = end

            for clustFts in clusters:
                if clustFts==[]:
                    continue
                else:
                    groupId = random.randint(0,maxInt)
                    seqIds = [f[seqIdx] for f in clustFts]
                    for ft in clustFts:
                        start, end, seqid, score = ft
                        if len(clustFts)!=1:
                            fusedSeqs.add(seqid) # keep a list of seqids that are fused with anything else
                        groupedFeatures.append( [chrom, start, end, seqid, score, seqIds, groupId] )
            
        # search and eliminate single features of seqids that were fused to another sequence elsewhere
        # XX due to feature reallocation this should not be necessary anymore 
        filtGrpFts = []
        for ft in groupedFeatures:
            chrom, start, end, seqid, score, seqIds, groupId = ft
            if len(seqIds)==1 and seqIds[0] in fusedSeqs:
                continue
            else:
                filtGrpFts.append(ft)
        return filtGrpFts

    #def predictedFeatures(self, maxDist):
        #""" convenience method: return clustered-bed features from filtered features """
        #bestGenome = self.bestGenome_mostSequences()
        #fts = self.chain(bestGenome, maxDist)
        #return bestGenome, fts

    def byGenome(self):
        # XX this changed in 09/09: returns now genome as well
        """ return as dictionary: genome -> list of (seqid, genome, chrom, start, end, score) """
        hits=self
        byGenome = {}
        for seqId, seqHits in hits.iteritems():
            for h in seqHits:
                genome, chrom, start, end, score, percId = h
                byGenome.setdefault(genome, []).append((seqId, genome, chrom, start,end, score, percId))
        return byGenome

    #def onlyGenome(self, genome):
        #""" return as dictionary: seqid -> list of (genome, chrom, start, end, score) """
        ##for k,v in self.byGenome().iteritems():
            #print k,v
        #hitGenomes = self.byGenome()
        #if genome in hitGenomes:
            #hits = hitGenomes[genome]
        #else:
            #return None
        #bySeqId = {}
        #for tuple in hits:
            #chrom, start, end, seqId, score = tuple
            #bySeqId.setdefault(seqId, []).append((genome, chrom, start,end, score))
        #return bySeqId

    def byChrom(self, genome):
        """ return as dictionary: chrom -> list of (start, end, seqid, score) """
        byChrom = {}
        for seqId, seqHits in self.iteritems():
            for h in seqHits:
                hitGenome, chrom, start, end, score = h
                if hitGenome==genome:
                    byChrom.setdefault(chrom, []).append((start,end, seqId, score))
        return byChrom
    
    def asBed(self):
        lines = {}
        for genome, genomeHits in self.byGenome().iteritems():
            for h in genomeHits:
                chrom, start,end, seqId, score = h
                data = [chrom, str(start), str(end), "seq"+str(seqId),str(score)]
                lines.setdefault(genome, []).append("\t".join(data))
        return lines

    def asText(self):
        lines = []
        for genome, bedLines in self.asBed().iteritems():
            lines.append( "- Genome:"+idToGenome[genome])
            lines.extend(bedLines)
        return "\n".join(lines)

    def bestGenomes(self, pmcId, bestTargetPredictor):
        """ return dict with method -> bestGenome """

        def prettyDict(data):
            list = []
            for key, value in data:
                list.append(idToGenome[value]+":"+(str(key)))
            return ", ".join(list)

        bestGenomes = {}
        descMethodDict = {"Best hit"                             : bestTargetPredictor.bestGenome_bestHit,
                          "Maximum sum of scores"                : bestTargetPredictor.bestGenome_bestTotalSum,
                          "Maximum number of matching sequences" : bestTargetPredictor.bestGenome_mostSequences
                }
        if len(self)!=0:
            for desc, bestGenomeFunction in descMethodDict.iteritems():
                # determine score according to method
                genomeScores = bestGenomeFunction(pmcId, self)
                # use only genomes with highest score
                bestScoreGenomes = pmcArticle.bestScoreElements(genomeScores)
                # weight these best genomes
                bestGenome, bestWeightGenomes, genomeWeights = bestTargetPredictor.weightGenomes(bestScoreGenomes)
                # and keep the details of the process
                bestGenomes[desc] = (bestGenome, genomeScores, bestWeightGenomes, genomeWeights)
        return bestGenomes

class PmcArticle:
    """ an object representing a PMC article """
    def asFasta(self):
        lines = []
        for seqId, (filename, seq,raw) in self.seqs.iteritems():
            lines.append( ">PMC"+str(self.pmcId)+"|"+str(seqId)+" fileName="+filename+"\n"+seq)
        return "\n".join(lines)

    def readFromMysql(self, db, numId, withSeqs=True):
        """ read from mysql database, return None if not found in Db """
        #pmcId = pmcId.strip("pmcA")
        self.numId = numId

        if withSeqs:
            self.seqs = {}
            recs = util.sql(db, "select fileName, seqId, sequence, context from sequences where numId='%s'" % str(numId))
            #if len(recs)==0:
                #return None
            for r in recs:
                filename, seqId, seq, raw = r
                self.seqs[seqId]=(filename, seq, raw)

        row = util.sql(db, "select pmid, pmcId from refs where numId='%s';" % str(numId))
        self.pmid, self.pmcId = row[0]

        recs = util.sql(db, "select * from bestGenomes where numId='%s';" % str(numId))
        if len(recs)>0:
            self.bestGenome = recs[0][1]
        else:
            self.bestGenomes= None

        recs = util.sql(db, "select genomeId, chrom, start, end, numId, seqid, score, groupId from features where numId='%s'; " % str(numId))
        self.chainHits = recs

        #recs = util.sql(db, "select * from blastHits where numId='%s'" % str(pmcId))
        #self.rawHits = IndexedBlastHits(recs)
        #recs = util.sql(db, "select * from filtHits where numId='%s'" % str(pmcId))
        #self.filtHits = IndexedBlastHits(recs)
        #recs = util.sql(db, "select * from blastHitsGenes where numId='%s'" % str(pmcId))
        #self.filtGeneHits = IndexedBlastHits(recs)

        return True

