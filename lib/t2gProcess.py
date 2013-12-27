# routines to process blat/blast matches into taxons and features thereon

import logging, os, random
import maxbio, maxTables, namedtuple, tabfile, util

logger = logging.getLogger("t2gProcess")

# --- MAIN FUNCTION IN THIS MODULE ---

def filterProtein(compBlastFile, proteinTableFilename):
    fieldParser = maxTables.TableParser(fileObj=compBlastFile, fileType="blastConvert") # create field name converter
    blockReader = maxTables.BlockReader(compBlastFile, 0) # split file on first field (=pmcId)

    protTable = open(proteinTableFilename, "w")
    headers = ["#documentId", "sequenceId", "best match in PDB", "BLAST score\n"]
    protTable.write("\t".join(headers))

    for (pmcId, block) in blockReader.readNext():
        pmcId = int(pmcId)
        recs = fieldParser.parseBlock(block)
        bestHits = removeNonOptHits(pmcId, recs)
        for bh in bestHits:
            data = [str(bh.pmcId), str(bh.seqId), bh.chrom, str(bh.score)]
            protTable.write("\t".join(data)+"\n")
        

def filter(compBlastFile, genomeWeights, helperTaxons, bestTaxons, geneWeights, genomeFeatureFile, geneFeatureFile, bestGenomeFile, bestGeneFile, maxDist):
    """ combined genome/gene processing """
    fieldParser = maxTables.TableParser(fileObj=compBlastFile, fileType="blastConvert") # create field name converter
    blockReader = maxTables.BlockReader(compBlastFile, 0) # split file on first field (=pmcId)

    bestGenomePredictor = mostMatchesWeightedPredictor(genomeWeights, helperTaxons, bestTaxons)

    docCount = 0
    for (pmcId, block) in blockReader.readNext():
        docCount +=1
        if docCount % 1000 ==0:
            logging.info("Processed %d documents" % docCount)
        pmcId = int(pmcId)
        recs = fieldParser.parseBlock(block)
        # keep only best hits
        bestHits = removeNonOptHits(pmcId, recs)
        if len(bestHits)==0: # everything matched univec
            logger.debug("docId %s: everything seems to match univec" % str(pmcId))
            continue

        bestGenomeInfo = bestGenomePredictor.getGenomeId(pmcId, bestHits) 
        if bestGenomeInfo==None:
            logger.debug("no best genome found in best hits, trying all hits")
            bestGenomeInfo = bestGenomePredictor.getGenomeId(pmcId, recs)  # shouldn't this be bestHits instead of recs?
            bestHits = recs
            if bestGenomeInfo==None:
                logger.debug("no best genome found, giving up")
                continue

        # keep only hits from best genome and separate them into genes/genomes
        genomeHits, geneHits = filterSplitHits(bestHits, bestGenomeInfo.bestGenome)
        
        # chain and disambiguate chains
        genomeFeatures  = disambiguateChains(chainHitsGenomes(genomeHits, maxDist=maxDist))
        rawGeneFeatures = disambiguateChains(chainHitsGenes(geneHits))

        # disambiguate genes
        bestGenes, geneFeatures = filterGeneFeatures(rawGeneFeatures, geneWeights)

        # write to output files, only if we have some supporting evidence
        if len(genomeFeatures)>0 or len(bestGenes)>0:
            writeTsvFeatures(genomeFeatureFile, pmcId, genomeFeatures)
            writeTsvFeatures(geneFeatureFile, pmcId, geneFeatures)
            writeTuples(bestGeneFile, pmcId, bestGenes)
            writeTuples(bestGenomeFile, pmcId, [bestGenomeInfo.bestGenome])

# ---- end main functions

#
# helper functions 
#

def filterGeneFeatures(features, geneWeights):
        """ some additional filters for genes: keep only genes with >=2 matches """

        # make dict: geneId -> matching seqIds
        seqIdGenes = {}
        for f in features:
            geneId = f.chainId
            seqIdGenes.setdefault(geneId, set()).add(f.seqId)

        logger.debug("Identified %d genes in results: %s" % (len(seqIdGenes), str(seqIdGenes)))

        # keep only geneIds with >=2 seqIds
        goodGenes = set()
        for geneId, seqIds in seqIdGenes.iteritems():
            if len(seqIds)>=2:
                goodGenes.add(geneId)

        logger.debug("Identified %d good genes (%s) with >=2 matches" % (len(goodGenes), str(goodGenes)))

        # keep only features on these genes
        filtFeatures = [f for f in features if f.chainId in goodGenes]
        logger.debug("Retained %d features on >=2 matching genes" % (len(filtFeatures)))

        # for all genes with identical sequence ids: take only the one with the best weight
        # [chainedFeature(genomeId=3702, chrom='AT1G67220.1-TAIR', start=1875, end=1896, seqId=2, score=44.0, seqIds=[2, 3, 4, 5, 6, 7], chainId='AT1G67220-TAIR-G'), chainedFeature(genomeId=3702, chrom='AT1G67220.1-TAIR', start=2772, end=2795, seqId=3, score=48.0, seqIds=[2, 3, 4, 5, 6, 7], chainId='AT1G67220-TAIR-G')]

        # index by seqIds: seqIdString -> genes
        seqIdGenes = {}
        for feat in filtFeatures:
            seqIdString = str(feat.seqIds)
            gene = feat.chainId
            seqIdGenes.setdefault(seqIdString, set()).add(gene)

        # for each seqIdString, take only the highest weighted gene
        bestGenes = set()
        for seqIdString, genes in seqIdGenes.iteritems():
            seqIdBestGenes, weightInfo = weightTargets(geneWeights, genes, strip=True)
            logger.debug("SeqIds %s: Ranking of genes got %s" % (seqIdString, str(weightInfo)))
            logger.debug("Best genes are %s" % seqIdBestGenes)
            bestGenes.update(seqIdBestGenes)

        return bestGenes, filtFeatures

def writeTuples(file, id, dataList):
    """ appends lines with id<tab>data<newline> to file """
    for number in dataList:
        row = [str(id)]
        row.append(str(number))
        line = "\t".join(row)+"\n"
        file.write(line)

#def chainDisambiguateHits(bestGenomeHits, maxDist, areGenes):
    # fuse neighboring hits into "chains" of hits, assign chainIds 
    #if isGenes:
        #chains = chainHitsGenes(bestGenomeHits)
    #else:
        #chains = chainHitsGenome(bestGenomeHits, maxDist) 
    # keep non-unique hits only in the chain with the most hits
    #features = disambiguateChains(chains) 
    #return features

def filterSplitHits(hits, bestGenomeId):
    """ keep only hits from best genome and split into genome/gene hits """
    # keep only hits on best genome
    bestGenomeHits = [hit for hit in hits if hit.genomeId==bestGenomeId] 
    logger.debug("%d hits on best genome" % (len(bestGenomeHits)))

    # separate bestGenomeHits into genome and gene bestGenomeHits
    geneHits = []
    genomeHits = []
    for h in bestGenomeHits:
        if "|" in h.chrom:
            geneHits.append(h)
        else:
            genomeHits.append(h)
    logger.debug("%d hits on best genome, %d hits on genes of best genome" % (len(genomeHits), len(geneHits)))
    return genomeHits, geneHits

def removeNonOptHits(pmcId, recs):
    logger.debug("START, document id %d" % pmcId)
    logger.debug("%d raw hits" % len(recs))
    bestHits = maxbio.bestTuples(recs, "seqId", "score") # keep only best matches
    logger.log(5, "%d best hits: %s" % (len(bestHits), bestHits))
    bestHits = removeUnivec(bestHits) # remove seqIds that match univec
    return bestHits

def writeTsvFeatures(file, pmcId, features):
    """ write tuples to tabsep file in format (genomeId, chrom, start, end, pmcId, seqId, score, chainId) """
    for f in features:
        data = [f.genomeId, f.chrom, f.start, f.end, pmcId, f.seqId, f.score, f.chainId]
        data = [str(e) for e in data]
        file.write("\t".join(data)+"\n")

def removeUnivec(bestHits):
    """ remove all hits with a genomeid==0 """
    noiseSeqIds = set([bh.seqId for bh in bestHits if bh.genomeId==0])
    result = [bh for bh in bestHits if bh.seqId not in noiseSeqIds]
    logger.debug("%d non-univec hits" % (len(bestHits)))
    if len(result)==0:
        logger.debug("No hits left, skipping this document")
    return result

ChainedFeature = namedtuple.namedtuple("chainedFeature", "genomeId, chrom, start, end, seqId, score, seqIds, chainId")

def chainHitsGenes(hits):
    """ index by geneId, convert to "chainedFeatures" (the chain is the just
    the gene here) and set the seqIds attribute to all seqids that match a
    particular gene """ 
    
    # index by gene
    logger.debug("%d unchained gene hits" % len(hits))
    geneFeatures = {}
    for f in hits:
        gene = f.chrom.split("|")[1]
        geneFeatures.setdefault(gene, []).append( f)

    # convert to chained features: chrom=transcript, chainId=gene
    chainedFeatures = {}
    for gene, features in geneFeatures.iteritems():
        allSeqIds = list(set([f.seqId for f in features]))
        newfeatures = []
        for f in features:
            trans = f.chrom.split("|")[0]
            newfeatures.append(ChainedFeature(f.genomeId, trans, f.tStart, f.tEnd, f.seqId, f.score, allSeqIds, gene))
        chainedFeatures[gene] = newfeatures

    logger.debug("converted to pseudo-chained features on %d genes" % len(chainedFeatures))
    return chainedFeatures


def chainHitsGenomes(hits, maxDist):
    """ fuse features into chains that are on same chromosome
    and closer than maxDist bp and tag all of them with a unique chainId.
    """
    logger.debug("%d unchained genome hits" % len(hits))

    maxInt = 2**32

    startIdx, endIdx, seqIdx, scoreIdx = 0, 1, 2, 3

    chromFts = maxbio.indexByField(hits, "chrom")

    chainedFeatures = {}
    # try to chain all features per chromosome
    for chrom, features in chromFts.iteritems():
        #print "new chrom", chrom
        #logger.debug("chaining: new chrom")
        features = maxbio.sortList(features, "tStart", reverse=False)
        chains = []
        lastEnd = None

        # if close enough: add to the last chain of features
        # if too far: break chain and start new chain
        for f in features:

            #print "feature", f
            #print "chains:", chains
            #if len(chains)>0:
                #print "lastchain", chains[-1]
            #print "pos", f.tStart, lastEnd 
            #if lastEnd!=None:
                #print "diff", f.tStart - lastEnd

            if lastEnd==None or abs(f.tStart - lastEnd) > maxDist: 
                #print "append empty one"
                chains.append([])
            chains[-1].append(f)
            #print "new last end", f.tEnd
            lastEnd = f.tEnd

        # mark all features with a unique chains id 
        # & store chainId as an attribute in all features of a chain
        for chainHits in chains:
            if chainHits==[]:
                continue
            else:
                chainId = random.randint(0,maxInt) # chainId is a random integer
                seqIds = [f.seqId for f in chainHits]
                chain = []
                for f in chainHits:
                    chain.append( ChainedFeature(f.genomeId, f.chrom, f.tStart, f.tEnd, f.seqId, f.score, seqIds, chainId) )
                chainedFeatures[chainId] = chain
    logger.debug("%d chained hits" % len(chainedFeatures))
    return chainedFeatures
        
def weightTargets(weights, entities, strip=False):
    """ given a list of entities and weights, return the ones with the highest weights""" 

    genomeWeights = []
    EntityWeight = namedtuple.namedtuple("entityWeight", "entityId, weight")
    for e in entities:
        score = weights.get(e,0) # if not in file, use score 0
        genomeWeights.append( EntityWeight(e,score) )
    bestElements = maxbio.bestScoreElements (genomeWeights, "weight")

    if strip:
        bestElements = [el.entityId for el in bestElements] # strip off weight information
    return bestElements, genomeWeights

def disambiguateChains(chains): 
    """ given a dict chainId -> list of members, returned a filtered list where members are 
    mapped ONLY to those chains with the most members. """

    # create weights for chains: number of members
    chainWeights = {}
    for chainId, hits in chains.iteritems():
        chainWeights[chainId] = len(hits)
    # index: seqId -> chains
    chainedHits = maxbio.flattenValues(chains)
    seqIdChains = maxbio.indexByField(chainedHits, "seqId")

    # for each seqId, determine chainIds(s) with highest weight
    bestChainIds = {}
    for seqId, hits in seqIdChains.iteritems():
        chainIds = set([h.chainId for h in hits])
        bestChains, chainWeightDetails = weightTargets(chainWeights, chainIds, strip=True)
        bestChainIds[seqId] = bestChains

    # keep only hits that contain bestChainIds
    bestFeats = []
    for hit in chainedHits:
        if hit.chainId in bestChainIds[hit.seqId]:
            bestFeats.append(hit)
    logger.debug("%d hits after reallocation between chains" % (len(bestFeats)))
    return bestFeats

# --- two classes that both return the best genome, either based on the hits or on a textfile (from a previous run)

# return-type of getGenome(...)
BestGenomeInfo = namedtuple.namedtuple("BestGenomeInfo", "filterTaxons, bestGenome, genomeScores, bestScoreWeights")

# interface
class bestTargetPredictor:
    """ determine the best matching genome/gene """
    def getGenomeId(self, pmcId, bestHits):
        pass
    
# textfile reader
class textfilePredictor(bestTargetPredictor):
    """ return best genome as obtained from <documentId>tab<genomeId> textfile """
    def __init__(self, bestGenomeFile):
        self.bestGenomes = tabfile.slurpdict(bestGenomeFile, asInt=True, keyAsInt=True)

    def getGenomeId(self, pmcId, bestHits):
        """ given the document id and hits indexed by seqId, return the best genome, info about how disambiguation selected it is written to logger at debug level"""

        genomeId = self.bestGenomes.get(pmcId, 0)
        if genomeId==0:
            logger.debug("document ID %s, no best genome found in textfile" % pmcId)

        bestGenomeInfo = BestGenomeInfo(bestGenome=genomeId, genomeScores={}, bestScoreWeights={})
        return bestGenomeInfo

# best matches predictor
class mostMatchesWeightedPredictor(bestTargetPredictor):
    """ given a list of hits, return he genome with the most matching sequences and the best weight """

    def __init__(self, weights, helperTaxons, bestTaxons):
        # list of best taxons determined by 
        # select taxId, genome, count(*) as count from entrezGene JOIN genomes ON genomes.genomeId=entrezGene.taxId group by taxId order by count;
        # on the text2genome database
        # They are: Arabidopsis thaliana Sorghum bicolor Drosophila yakuba Danio rerio Caenorhabditis elegans Populus trichocarpa Oryza sativa Rattus norvegicus Saccharomyces cerevisiae Drosophila melanogaster Homo sapiens Mus musculus
        self.weights = weights
        self.helperTaxons = helperTaxons
        self.bestTaxons = bestTaxons

    @staticmethod
    def genomeSeqScore(indexedHits):
        """ count the number of sequences that match to each genome 
            and return as list of (taxonId, score) """

        genomeScores = {}
        # for each genome, count how many seqIds match
        for seqId, seqHits in indexedHits.iteritems():
            genomes = set([h.genomeId for h in seqHits])
            for g in genomes:
                genomeScores.setdefault(g, 0)
                genomeScores[g]+=1

        genomeScores = genomeScores.items() # convert to list (genome, score)
        GenomeScore = namedtuple.namedtuple("GenomeScore", "taxId, score")
        genomeScores = [GenomeScore(*g) for g in genomeScores] # created named tuples

        maxbio.sortList(genomeScores, "score") # sort, only to make the debug message easier to read
        logger.debug("genomeScores %s" % str(genomeScores))
        return genomeScores

    def getGenomeId(self, pmcId, bestHits):
        """ given the document id and hits indexed by seqId, return the best genome, info about how disambiguation selected is written to logger at debug level"""

        logger.debug("start determining best genome: %d raw hits" % (len(bestHits)))
        # if we have helper information: keep only matches from the helper taxons
        if pmcId in self.helperTaxons:
            filterTaxons = self.helperTaxons[pmcId]
            logger.debug("restricting matches to helper taxons: %s" % str(filterTaxons))
        # otherwise use only top10 genomes from entrez gene
        else:
            filterTaxons = self.bestTaxons
            logger.debug("No helper taxons, using only blast matches on best taxons: %s" % str(filterTaxons))
        logger.log(5, "hits before filtering on filterTaxons:"+str(bestHits))
        bestHits = [x for x in bestHits if x.genomeId in filterTaxons]
        logger.log(5, "hits after filtering on filterTaxons:"+str(bestHits))
        if len(bestHits)==0:
            logger.debug("Zero hits left, exiting")
            return None
         
        indexedHits = maxbio.indexByField(bestHits, "seqId")
        logger.debug("determining best genome: %d hits of %d sequences" % (len(bestHits), len(indexedHits)))

        # select the genomes with the most matching sequences
        genomeScores = self.genomeSeqScore(indexedHits)
        bestScoreGenomes = [s.taxId for s in maxbio.bestScoreElements(genomeScores, 1)]
        logger.debug("document %s, bestScoringGenomes %s" % (str(pmcId), str(bestScoreGenomes)))
        
        if len(bestScoreGenomes)!=1:
            # dis-ambiguate genomes
            # choose the genome with the highest taxon weight
            logger.debug("Disambiguation needed")
            bestWeightGenomes, bestScoreWeights = weightTargets(self.weights, bestScoreGenomes)
            if len(bestWeightGenomes)!=1:
                #random.shuffle(bestWeightGenomes)
                bestWeightGenomes.pop()
                logger.debug( "warning: document %s, had to select best genome randomly" % str(pmcId))
            bestWeightGenome = bestWeightGenomes[0].entityId
            logger.debug( "document %s: best disambiguated genomes %s" % (str(pmcId), str(bestWeightGenomes)))
        else:
            logger.debug("No disambiguation necessary")
            bestWeightGenome = bestScoreGenomes.pop()
            bestScoreWeights = []
        bestGenomeInfo = BestGenomeInfo(bestGenome=bestWeightGenome, filterTaxons=filterTaxons, genomeScores=genomeScores, bestScoreWeights=bestScoreWeights)
        return bestGenomeInfo

def trackBestReferences(pmcId, refGenomes, bestGenomeInfo):
    """ write to logger where the best genome is located in bestGenomeInfo """
    if not refGenomes:
        return

    # several genomes with topranked score  -> weighting was used
    weighting = False
    if len(maxbio.bestScoreElements(bestGenomeInfo.genomeScores, scoreField="score"))>1:
        weighting=True

    # prep scores
    scores = list(set([rec.score for rec in bestGenomeInfo.genomeScores]))
    scores.sort(reverse=True)
    # prep weights
    weights = list(set([entity.weight for entity in bestGenomeInfo.bestScoreWeights]))
    weights.sort(reverse=True)

    # identify ranks of ref genomes scores and ref genome weight in all scores/weights
    refRanks = []
    for refGenomeId in refGenomes:
        # find score of ref genome in scorse
        refGenomeScores = [rec.score for rec in bestGenomeInfo.genomeScores if rec.taxId==refGenomeId]
        if len(refGenomeScores)==0:
            refGenomeRank = -1 # not found
        else:
            refGenomeScore = refGenomeScores.pop()
            refGenomeRank = scores.index(refGenomeScore)
        refRanks.append(str(refGenomeRank))

    if weighting==True:
        refWeightRanks = []
        for refGenomeId in refGenomes:
            # find weight of ref genome
            refGenomeWeights = [rec.weight for rec in bestGenomeInfo.bestScoreWeights if rec.entityId==refGenomeId]
            if len(refGenomeWeights)==0:
                refGenomeRank=-1
            else:
                refGenomeWeight=refGenomeWeights.pop()
                refGenomeRank = weights.index(refGenomeWeight)
            refWeightRanks.append(str(refGenomeRank))
        refWeightString = ",".join(refWeightRanks) 
    else:
        refWeightString = "NoWeighting"

    data = [pmcId, ",".join([str(x) for x in refGenomes]), str(bestGenomeInfo.bestGenome), str(int(weighting)), ",".join(refRanks), refWeightString]
    print "\t".join(data)


if __name__ == "__main__":
    import doctest
    doctest.testmod()


