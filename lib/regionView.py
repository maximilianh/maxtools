# functions to export bed regions to html and to scan them for a given set of annotations

import sys
import html
import tabfile
import util
import math
import util
import types

# little helpers
def getGenes(b):
    """ get the geneLeft/geneRight of a bed entry, can use either block-annotated regions or old-style bed compatible annotation """
    if "block" in b.__dict__:
        #gene1 = b.block.geneLeft.split("_")[0]
        #gene2 = b.block.geneRight.split("_")[0]
        #genes, intronic = b.block.getGenes(flankType, maxDist)
        genes = b.flankGenes
    else:
        genes = []
        genes.append( b.name.split("|")[1].split("_")[0] )
        genes.append( b.name.split("|")[2].split("_")[0] )
    return genes

class RegionView:
    """ input bed-regions, output html for them"""
    def __init__(self, beds, annotations, htmlFh, geneScores, db, track, quiet, manualAnnotations, ishUrl, geneUrl):
        """ beds are beds with (optional) attributes (geneLeft, geneRight, alternative: in name-field of bed, split by |),
            annotations is a string with keyword=>filename (format: ghost=test.txt,c=d,...) , 
            where keyword is one of jgi1,aniseed, etc..."""
        """ geneScores is dict geneId -> length of conserved sequences around """
        self.quiet = quiet
        self.beds = beds
        self.annotations = {}
        self.of = htmlFh
        self.db = db
        self.track = track
        self.geneScores = geneScores
        self.ishUrl = ishUrl
        self.geneUrl = geneUrl

        allowedKeys = ["hugo", "human", "jgi1", "aniishsobral", "ghostish", "zfinish", "aniishfabrice"]
        if annotations!=None:
            self.debug("htmlOutput reading annotation files: ")
            annotations = dict([(f.split("=")[0], f.split("=")[1]) for f in annotations.split(",")])

            for key,fname in annotations.iteritems():
                key=key.lower()
                self.debug("%s " % fname)
                if key not in allowedKeys:
                    self.debug("error: key %s is not a valid keyword\n" % key)
                    sys.exit(1)
                d = tabfile.slurpdictlist(fname)
                self.annotations[key]=d
        self.debug("\n")
        # inverse annotation:
        annots = {}
        for annot, genes in manualAnnotations.iteritems():
            for g in genes:
                annots.setdefault(g, []).append(annot)
        self.annotations["activeAnnot"]=annots

    def debug(self, msg):
        if not self.quiet:
            sys.stderr.write(msg)

    def _write(self, text):
        self.of.write(text+"\n")

    def _writeHeader(self,title="Matches"):
        self.of.write("""
        <html>
        <head>
        <title>"""+title+"""</title>
        """)
        self._writeStylesheet()
        self._write("</head>")
        self._write("<body>")

    def _writeFooter(self,title="Matches"):
        self.of.write("""
        </body>
        </html>
        """)

    def _indexByGene(self,beds):
        idx = {}
        for b in beds:
            genes = getGenes(b)
            for gene in genes:
                idx.setdefault(gene, []).append(b)
        return idx

    def _sortGenesByMatchCount(self,index):
        genes = index.keys()
        genes.sort(key=lambda x: len(index[x])/self.geneScores[x], reverse=True)
        return genes

    def _writeStylesheet(self):
        self._write("""
        <style type="text/css">
        table {
                font: 10px Verdana, Arial, Helvetica, sans-serif;
                border-collapse: collapse;
                table-layout: fixed; 
              }
         td {
                padding: 2px;
                vertical-align: top;
                text-align: left;
            } 


        * { font-family:Arial,sans-serif; font-weight:normal; }
        .geneHead { width=200px; font-size:12px; font-weight:bold; background-color: #DDDDDD; border-top: 24px solid #FFFFFF; margin-top: 2em; }
        .geneGood    { background-color: green; font-weight: bold; }
        .geneNotGood { background-color: red; font-weight: bold; }
        .geneDetailsBold { background-color: grey; font-weight: bold; }
        .geneDetails { width=200px; font-size:11px; border-bottom: 1px solid #BBB; font-weight: bolder; }
        .regions     { font-size: 12px; }
        .regionHead  { font-size: 12px; font-weight:bold; background-color: #EEEEEE; }
        .exprText    { font-size:12px; font-weight: normal; text-align:left; padding-left: 10px; vertical-align: top; }
        </style>
        """)
 
    def _writeAnnotLinkListCell(self, keyword, gene, targetUrlFunc, style=None):
            
        if style!=None:
            self._write('<td class="%s">' % style)
        else:
            self._write("<td>")
        if keyword in self.annotations:
            links = []
            for sym in self.annotations[keyword].get(gene,["NotFound"]):
                if type(targetUrlFunc)==types.FunctionType:
                    targetUrl = targetUrlFunc(sym)
                else:
                    targetUrl = targetUrlFunc
                links.append(targetUrl)
            self._write(", ".join(links))
        else:
            if targetUrlFunc!=None:
                if type(targetUrlFunc)==types.FunctionType:
                    self._write(formatFunc(gene))
                else:
                    self._write('<A HREF="%s">%s</A>\n' % (targetUrlFunc, gene))
            else:
                self._write(gene)


        self._write("</td>")

    def _writeExprTextRow(self, title, orthoKey, annotKey, gene, stageFilt=None):
        # translate gene to orthologs
        if annotKey in self.annotations:
            if orthoKey in self.annotations:
                orthos = self.annotations[orthoKey].get(gene, [])
            else:
                orthos = [gene]
            for o in orthos:
                dict = {}
                dict = {}
                for text in self.annotations[annotKey].get(o, []):
                    #key = data
                    #clone,stage,text = data
                    #key = "|".join([clone,stage])
                    dict.setdefault(o, []).append(text)

                for id, texts in dict.iteritems():
                    #clone,stage = key.split("|")
                    text = ", ".join(texts)
                    self._write('<tr>')
                    self._write('<td class="exprText">%s</td>' % title)
                    title=""
                    self._write('<td class="exprText">%s</td>' % gene)
                    self._write('<td class="exprText">%s</td>' % o)
                    #self._write('<td>%s</td>' % stage)
                    self._write('<td class="exprText">%s</td>' % text)
                    self._write('</tr>')

    def _writeGene(self, gene, db, geneHasTargetAnnotation, geneHasNotTargetAnnotation):
            self._write('<tr class="geneHead">')
            self._write('<td class="geneHead">Gene Symbol</td>')
            self._write('<td class="geneHead">Gene Id</td>')
            #self._write('<td class="geneHead">Human Gene Id</td>')
            #self._write('<td class="geneHead">Aniseed Gene</td>')
            self._write('<td class="geneHead">InSitus</td>')
            self._write('</tr>')

            self._write('<tr class="geneDetails">')
            # col: hugo gene symbols + genecards links
            style = "geneDetailsBold"
            if gene in geneHasNotTargetAnnotation:
                style = "geneNotGood"
            elif gene in geneHasTargetAnnotation:
                style = "geneGood"
            self._writeAnnotLinkListCell("hugo", gene, html.geneCardsLink, style)
            # col: ensembl gene link
            #self._write("<td>%s</td>" % html.ensGeneLink(gene, db))
            self._write("<td><A HREF='%s'>%s</A></td>" % (self.geneUrl.replace("$GENE", gene), gene))
            # col: human gene, ensembl linked
            #self._writeAnnotLinkListCell("human", gene, html.ensGeneLink)
            # col: link to aniseed page
            #self._writeAnnotLinkListCell("jgi1", gene, html.aniGeneLink)
            # col: link to aniseed pictures
            #self._writeAnnotLinkListCell("jgi1", gene, html.aniISHLink)
            self._writeAnnotLinkListCell("jgi1", gene, self.ishUrl.replace("$GENE", gene))
            self._write("</tr>")

            # Rows for expression descriptions
            self._write("<tr>")
            self._write("<td colspan='5'>")
            self._write("<table>")
            #self._writeExprTextRow("Aniseed(Daniel)", "jgi1", "aniish", gene, ["Tadpole", "Mid tailbud", "Initial tailbud", "Late tailbud"])
            #self._writeExprTextRow("Active Annotation", "jgi1", "activeAnnot", gene)
            self._writeExprTextRow("Aniseed/Daniel", "jgi1", "aniishsobral", gene)
            self._writeExprTextRow("Aniseed/Fabrice", "jgi1", "aniishfabrice", gene)
            self._writeExprTextRow("Ghost", "jgi1", "ghostish", gene)
            self._writeExprTextRow("Zfin", "jgi1", "zfinish", gene)
            #self._writeExprTextRow("ZFIN", "dr", "zfin", gene)
            self._write("</table>")
            self._write("</td>")
            self._write("</tr>")

    def _writeRegions(self, regions, hgsid, baseUrl, db, track):
            self._write('<tr><td class="regionHead" colspan="3">Conserved regions around this gene that match your motifs:</td></tr>')
            for h in regions:
                self._write("<tr>")
                self._write("<td class='regions'>"+html.ucscGenomeLink(baseUrl, db, "%s:%d-%d" % (h.chrom, h.start, h.end), hgsid=hgsid)+"</td>")
                self._write("<td class='regions'>"+html.ucscMafLink("(alignment)", baseUrl, db, track, h.chrom, h.start, h.end)+"</td>")
                self._write("</tr>")

    def _writeStats(self, flankingGenes, geneCount, fgGenes, bgGenes):
        def writeCountPercentage(text, count, total):
                self._write("<li>")
                self._write(text)
                self._write(str(count))
                if total!=0:
                    self._write("(%.3f%% out of %d)" % (100 * float(count) / total, total))
                self._write("</li>")

        self._write("<ul>")

        counter = exprAnnotCounter(fgGenes, bgGenes, self.geneScores)
        genomeGenesWithAnnot, genomeGenesWithTAnnot = counter.genomeStats()
        hitStats = counter.hitStatsBeds(self.beds, noHyperg=True)
        annotGenomeGenesNotTAnnot = genomeGenesWithAnnot - genomeGenesWithTAnnot

        # output 
        self._write("<li>Genomic hits: %d</li>" % len(self.beds))
        #self._write("<li>Genes in genome: %d</li>" % geneCount)
        writeCountPercentage("Hits with at least one annotated flanking gene:",hitStats.regionsWithAnnot, len(self.beds))
        #writeCountPercentage("Number of annotated flanking genes:", hitStats.flankGenesWithAnnot, len(flankingGenes))
        writeCountPercentage("Hits with at least one flanking annotated gene expressed in target tissue:", hitStats.regionsWithTargetAnnot, hitStats.regionsWithAnnot)
        writeCountPercentage("Genes flanking the hits:",len(flankingGenes), geneCount)
        writeCountPercentage("Annotated genes in genome (background + foreground)" , genomeGenesWithAnnot, geneCount)
        writeCountPercentage("Annotated genes in genome expressed in target tissues:", genomeGenesWithTAnnot, genomeGenesWithAnnot)
        #writeCountPercentage("Total number of annotated genes in genome not expressed in target regions (%s):", annotGenomeGenesNotTAnnot, genomeGenesWithAnnot)
        self._write("<li>Annotated genes flanking a hit: %d</li>" % len(hitStats.flankingAnnotatedGenes))
        self._write("<li>Annotated genes flanking a hit AND expressed in target region: %d</li>" % hitStats.TP)

        self._write("<li>Quality of hits:</li>")
        self._write("<ul>")
        self._write("<li>TP=%.2f, TN=%.2f, FP=%.2f, FN=%.2f</li>" % (hitStats.TP, hitStats.TN, hitStats.FP, hitStats.FN))
        self._write("<li>Sensitivity=%f</li>" % (hitStats.Sens))
        self._write("<li>Specificity=%f</li>" % (hitStats.Spec))
        self._write("<li>Positive Predictive Value=%f</li>" % (hitStats.PPV))
        #self._write("<li>Performance Coefficient=%f</li>" % (hitStats.PC))
        self._write("<li>Pearson Correlation=%f</li>" % (hitStats.CC))
        #enrich = math.log10(hitStats.hypergParams)
        #self._write("<li>Enrichment score (Harbison et al Nature 2004, 1 - log(hyp-pValue)): %f</li>" % enrich

        self._write("</ul>")

        #self._write("<li>Hypergeometric Probability: There is a box with %(N)d balls of which %(m)d are white. One draws %(n)d balls. The probability to find %(k)d or more white balls is: %(pVal).10g</li>" % hitStats.hypergParams)
        self._write("<li>Binomial Probability: When repeating an experiment %(n)d times where the probability of success is %(p)f, the probability to obtain %(k)d or more 1s is: %(pVal).10g</li>" % hitStats.bnpParams)
        #self._write("<li>Poisson Probability: We throw a coin 0/1 %(n)d times. The probability to obtain a one is %(p)f. Lambda (=n*p) is %(lambda)g. The probability to obtain %(k)d  or more 1s is: %(pVal).10g</li>" % hitStats.poissParams)
        #self._write("<li>Corrected Binomial Probability: Around target genes we have on average %(avgConsTarget)d bp conserved, those around annotated genes are on average %(avgConsAnnot)d. So there are %(corrFactor)f times more CNS around targets than the average gene. Parameters for binomials are: n = %(n)d, p= %(p)f, k<=%(k)d. Resulting pValue is %(pVal).10g</li>" % hitStats.corr_bnpParams)

        #self._write("<li>Sensitivity (how well do the hits flank the positive genes): The hits flank %.2f percent of all genes with this annotation (%d of %d)</li>"% ((100 * float(len(genesWithTAnnotDict)) / genomeGenesWithTAnnot), len(genesWithTAnnotDict), genomeGenesWithTAnnot))
        #self._write("<li>Specificity (how well do the hits avoid the negative genes): %.2f percent of all annot genes without a target annotation (%d of %d) are not flanked by a hit</li>"% ((100 * float(len(annotNotTargetNotFlank)) / annotGenomeGenesNotTAnnot), len(annotNotTargetNotFlank), annotGenomeGenesNotTAnnot))
        self._write("</ul>")
        return counter


    def writeMatchesHtml(self, geneCount, fgGenes, bgGenes, motifMatchBeds, sortedGenes, geneToRegions, browserUrl, ucscDb, clade, organism):
        # upload track to ucsc servers
        data = str(motifMatchBeds)
        hgsid = html.ucscUpload(ucscDb, data, server=browserUrl, name="Motif matches", description="Motif matches", clade=clade, organism=organism)


        self._write("<h3>Overview of matches in genome:</h3>")
        statCounter = self._writeStats(sortedGenes, geneCount, fgGenes, bgGenes)
        self._write("<h3>Matches: Genes with flanking hits</h3>")
        self._write("Some notes the following list:<small><ul>")
        self._write("<li>grey = no expression data available</li><li>red = predicted but background annotation</li><li>green = predicted and foreground annotation</li>")
        self._write("<li>Clicking on genomic coordinates will link to a UCSC genome browser with the conserved motif matches highlighted with a custom track</li>")
        self._write('<li>These coordinates can be pasted into <a href="../genomePrimers/genomePrimers.cgi">Genome Primers</a> to create Gateway primers</li>')
        self._write("<li>Sorted by number of hits / length of flanking conserved bp</li>")
        self._write("<li>Gene annotation data from Zebrafish (from zfin.org, stage: 24hrs) annotation has been mapped with best-reciprocal blast hits but is displayed just for information</li>")
        self._write("<li>C. intestinalis: Ghost annotation data has been downloaded from Ghost and, GhostIds have been semi-manually mapped to JGI1-model IDs</li>")
        self._write("<li>C. intestinalis: Aniseed/Fabrice is the data as received from the Aniseed Group (Fabrice Daian) in 2008</li>")
        self._write("</ul></small>")

        self.writeHiddenForm(fgGenes, "Export hits (bed format)", type="bed")
        self.writeHiddenForm(fgGenes, "Export green genes", type="tp")
        self.writeHiddenForm(fgGenes, "Export red genes", type="fp")

        self._write("<table>")
        # write html
        bedLines = []
        for gene in sortedGenes:
            if gene=="Gap":
                continue
            hits = geneToRegions[gene]
            self._writeGene(gene, "Ciona_intestinalis", statCounter.geneHasTargetAnnot, statCounter.geneHasNotTargetAnnot)
            self._writeRegions(hits, hgsid, browserUrl, self.db, self.track)

        self._write("</table>")


    def writeStatsHtml(self, geneCount, fgAnnots, bgGenes, sortedGenes):
        """ output only general statistics about number of matches, quality, etc, not links to genome browser, no gene lists """

        self._write("<h3>Matches:</h3>")

        #flankAnnotated = len(set(sortedGenes).intersection(annotatedGenes))
        counter = exprAnnotCounter(set(['dummy']), bgGenes, self.geneScores)
        genomeGenesWithAnnot, genomeGenesWithTAnnot = counter.genomeStats()
        hitStats = counter.hitStatsBeds(self.beds, noHyperg=True)

        self._write("<ul>\n")
        self._write("   <li>Hits: %d</li>\n" % len(self.beds))
        self._write("   <li>Flanking genes: %d</li>\n" % len(sortedGenes))
        self._write("   <li>Flanking genes with any annotation: %d (out of %d total annotated genes)</li>\n" % (len(hitStats.flankingAnnotatedGenes), genomeGenesWithAnnot))
        self._write("</ul>\n")

        self._write("<bold>Please note: you have specified <i>several</i> tissues. As a result, a summary of p-Values per tissue will be shown and not the individual matches per motif-combination:</bold>")
        self._write("<h3>Quality of matches by annotated tissue:</h3>")
        self._write("<ul>\n")
        lines = []
        for annot, targetGenes in fgAnnots.iteritems():
            counter = exprAnnotCounter(targetGenes, bgGenes, self.geneScores)
            hitStats = counter.hitStatsBeds(self.beds, noHyperg=True)
            #lines.append( [hitStats.bnpParams["pVal"], ("<li>%s:<ul><li>%d genes with matches out of %d</li><li>binomial: %g</li><li>hyperg.: %g</li><li>Predictive Power: %g</li><li>Correlation: %g</li><li>corrected binomial prob: %g</li></ul></li>" % (annot, hitStats.TP, len(targetGenes), hitStats.hypergParams["pVal"], hitStats.bnpParams["pVal"], hitStats.PC, hitStats.CC, hitStats.corr_bnpParams["pVal"]))])
            lines.append( (hitStats.bnpParams["pVal"], ("<li>%s:<ul><li>%d genes with matches out of %d (share: %g)</li><li>Random background rate: %g</li><li>binomial: %g</li><li>Predictive Power: %g</li><li>Correlation: %g</li></ul></li>" % (annot, hitStats.TP, len(targetGenes), float(hitStats.bnpParams["k"])/hitStats.bnpParams["n"], hitStats.bnpParams["p"],hitStats.bnpParams["pVal"], hitStats.PC, hitStats.CC))))
        lines.sort(key=lambda (x,y): x)
        self._write("\n".join([y for x,y in lines]))
        self._write("</ul>\n")


    def writeHiddenForm(self, fgGenes, text, type):
        lines = set()
        # generate data to output on form
        if type=="bed":
            for b in self.beds:
                lines.add("\t".join([b.chrom, str(b.start), str(b.end), str(b.name.split("|")[0]), "0", "+"]))
        elif type=="tp" or type=="fp":
            for b in self.beds:
                genes = getGenes(b)
                for g in genes:
                    if type=="tp":
                        if g in fgGenes:
                            lines.add(g)
                    elif type=="fp":
                        if g not in fgGenes:
                            lines.add(g)

        # write form
        print '<form action="showRawData.cgi" method="post">'
        i = 0
        for l in lines:
            i+=1
            print ('<input type="hidden" name="data" value="%s">' % l)
        print '<input type="submit" value="%s">' % (text)
        print '</form>'

    def writeHtml(self, foregroundAnnotations, allAnnotatedGenes, geneCount, parameters, motifMatchBeds, browserUrl, ucscDb, clade, organism):

        self._writeHeader()

        # index matches by gene
        gene2Region = self._indexByGene(self.beds)
        sortedGenes = self._sortGenesByMatchCount(gene2Region)

        if len(foregroundAnnotations)==1:
            foregroundGenes = foregroundAnnotations.values()[0]
            self.writeMatchesHtml(geneCount, foregroundGenes, allAnnotatedGenes, motifMatchBeds, sortedGenes, gene2Region, browserUrl, ucscDb, clade, organism)
        else:
            self.writeStatsHtml(geneCount, foregroundAnnotations, allAnnotatedGenes, sortedGenes)

        self._write("<h3>Parameters used for this run:</h3>")
        self._write("<ul>")
        for name, val in parameters.iteritems():
            self._write("<li>")
            self._write("<b>"+name+":</b> "+str(val))
            self._write("</li>")
        self._write("</ul>")

        self._writeFooter()
        self.of.close()


# ==================================== 
class exprAnnotCounter:
    def __init__(self, fgGenes, annotGenes, geneScores):
        self.geneHasAnnot = annotGenes
        self.geneHasTargetAnnot = fgGenes
        self.geneHasNotTargetAnnot = annotGenes.difference(fgGenes)
        self.geneScores = geneScores

        if len(fgGenes)==0:
            sys.stdout.write("error: no annotation file specified and no target genes specified\n")
            sys.stderr.write("error: no annotation file specified and no target genes specified\n")
            sys.exit(1)

    def genomeStats(self):
        """ some general genome statistics, do not depend on the hits """
        genesWithAnnot = len(self.geneHasAnnot)
        genesWithTargetAnnot = len(self.geneHasTargetAnnot)
        return genesWithAnnot, genesWithTargetAnnot

    def _addPVal(self, stats, noHyperg=False):
        N    = len(self.geneHasAnnot)
        m    = len(self.geneHasTargetAnnot)
        n    = len(stats.flankingAnnotatedGenes)
        k    = stats.TP

        assert(N!=0)
        assert(m!=0)

        if noHyperg:
            pVal_hgm=1.0
        else:
            pVal_hgm = 1.0 - util.hypergProbSum(k, N, m, n)
            stats.hypergPval    = float(pVal_hgm)
            stats.hypergParams  = {'N': N, 'm' : m, 'n' : n, 'k' : k, 'pVal' : pVal_hgm }

        if N!=0:
            p = float(m)/N
        else:
            p = 0.0

        # binom Probability
        # we need to lower k by 0.5 as we calculate the p that x > k 
        # in the case of a 100% correct prediction, p will be zero and the log undefined
        # we avoid this by calculating p that x > k - 0.5
        # and we avoid the decision to which class the point mass at k belongs (huh?)
        pVal_bp = util.binProbGt(k-0.5, size=n, prob=p) 
        stats.bnpPval = pVal_bp
        if pVal_bp!=0.0:
            stats.bnpScore = -math.log10(pVal_bp)
        else:
            stats.bnpScore = 1000000000 # super-duper score, should never be obtained by any prediction

        stats.bnpParams  = {'n': n, 'k' : k, 'pVal' : pVal_bp, 'p' : p, 'score' : stats.bnpScore}

        # poisson probability 
        #sys.stderr.write("n, p, k: %d, %f, %d\n" % (n, p, k))
        #pVal_poiss = 1.0 - util.poissProbSum(n, p, k) 
        #sys.stderr.write("old P-Val: %f\n" % pVal_poiss)
        #pVal_poiss = 1.0 - util.poissProbSum_scipy(n, p, k) 
        pVal_poiss = 1.0 - util.poissProbSum(n, p, k) 
        #sys.stderr.write("new P-Val: %f\n" % pVal_poiss)
        stats.poissParams  = {'lambda' : n*p, 'n': n, 'k' : k, 'pVal' : pVal_poiss, 'p' : p}

        # corrected binom. probab., using relation target CNS len / all CNS len as p
        # using priors of 0.01 to avoid 0 probabilities
        geneScores = self.geneScores
        targetScore       = sum([geneScores.get(g,0.1) for g in self.geneHasTargetAnnot]) 
        annotScore        = sum([geneScores.get(g,0.1) for g in self.geneHasAnnot]) 
        flankAnnotScore   = sum([geneScores.get(g,0.1) for g in stats.flankingAnnotatedGenes]) 
        flankTargetScore  = sum([geneScores.get(g,0.1) for g in stats.flankingTargetGenes]) 
        avg_All_Score = float(annotScore)/  N
        avg_Trg_Score = float(targetScore)/ m

        corrFactor = (avg_Trg_Score / (avg_All_Score+1))
        corr_p =  corrFactor * p

        if k>0:
            corr_pVal_bp = util.binProbGt(k, size=n, prob=corr_p) 
        else:
            corr_pVal_bp = "No matches"

        stats.corr_bnpPval    = corr_pVal_bp
        stats.corr_bnpParams  = {'consTarget': targetScore, 'consAnnot' : annotScore, 'consFlankAnnot' : flankAnnotScore, 'consFlankTarget' : flankTargetScore, 'avgConsTarget' : avg_Trg_Score, 'avgConsAnnot' : avg_All_Score,'n': n, 'k' : k, 'pVal' : corr_pVal_bp, 'p' : corr_p, 'corrFactor' : corrFactor}

    def hitStatsGenes(self,predictedGenes, pVal=True, noHyperg=False):
        """ get stats given a set of predicted genes """
        stats = util.hitStats(self.geneHasAnnot, predictedGenes, self.geneHasTargetAnnot)
        stats.flankingAnnotatedGenes           = predictedGenes.intersection(self.geneHasAnnot)
        stats.flankingTargetGenes              = predictedGenes.intersection(self.geneHasTargetAnnot)
        self._addPVal(stats, noHyperg)
        return stats

    def hitStatsBeds(self,beds, pVal=True, noHyperg=False):
        """ get stats given a list of bed features linked to wordBlocks, linked to genes """
        regionsWithAnnot = 0
        regionsWithTargetAnnot = 0
        predictedGenes = set()

        for b in beds:
            genes = getGenes(b)
            for gene in genes:
                predictedGenes.add(gene)
                if gene in self.geneHasTargetAnnot:
                    regionsWithTargetAnnot +=1          # we assume that if gene hasTargetAnnot then geneHasAnnot
                if gene in self.geneHasAnnot: 
                    regionsWithAnnot +=1
                    break                                # we only want to count each region once
        
        stats = self.hitStatsGenes(predictedGenes, pVal, noHyperg=noHyperg)
        stats.regionsWithAnnot = regionsWithAnnot
        stats.regionsWithTargetAnnot = regionsWithTargetAnnot

        return stats

