import collections, subprocess
from os.path import basename
import json

# convert icgc to knowledger submission format

def writeMeta(projCode, projMetaDict):
    " write project meta data json file "
    subMetaFname = "ledgerSubmit/%s.json" % projCode
    print "writing meta to %s" % subMetaFname
    json.dump(projMetaDict, open(subMetaFname, "w"), sort_keys=True, indent=4, separators=(',', ': '))

validChroms = ["X", "Y", "1", "2","3","4","5", "6","7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "MT"]

# the fields from the simple somatic mutation files that we keep in the sample metadata
expFields = tuple(["gene_build_version", "platform","experimental_protocol", "sequencing_strategy", "base_calling_algorithm", "alignment_algorithm", "variation_calling_algorithm", "other_analysis_algorithm", "seq_coverage", "raw_data_repository", "raw_data_accession", "initial_data_release_date"])

vcfHead = """##fileformat=VCFv4.1
##comment=ICGC open access Simple Somatic Mutations (SSM) in VCF format, converted by icgcToLedger.py
##reference=GRCh37
##source=ICGC20
##INFO=<ID=mutated_from_allele,Number=1,Type=String,Description="original allele in sequenced tumor can be different from ref">
##INFO=<ID=TC,Number=1,Type=String,Description="total read count">
##INFO=<ID=MC,Number=1,Type=String,Description="mutant allele read count">
##INFO=<ID=GA,Number=1,Type=String,Description="gene affected">
##INFO=<ID=TA,Number=1,Type=String,Description="transcript affected">
##INFO=<ID=VS,Number=1,Type=String,Description="validation status">
##INFO=<ID=VP,Number=1,Type=String,Description="validation platform">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
"""

import glob, re, collections, gzip

def nextRow(inFile, encoding="utf", fieldSep="\t"):
    """ 
    simplified version of iterTsvRows, return rows from tab-sep file
    """
    fh = inFile

    line1 = fh.readline()
    line1 = line1.strip("\n").strip("#")
    headers = line1.split("\t")
    #headers = [re.sub("[^a-zA-Z0-9_]","_", h) for h in headers]

    Record = collections.namedtuple('tsvRec', headers)
    for line in fh:
        line = line.strip("\n")
        fields = line.split(fieldSep)
        if encoding!=None:
            fields = [f.decode(encoding) for f in fields]
        try:
            rec = Record(*fields)
        except Exception, msg:
            logging.error("Exception occured while parsing line, %s" % msg)
            logging.error("Filename %s" % fh.name)
            logging.error("Line was: %s" % line)
            logging.error("Does number of fields match headers?")
            logging.error("Headers are: %s" % headers)
            raise Exception("wrong field count in line %s" % line)
        # convert fields to correct data type
        yield rec

print "reading proj codes"
projectInfo = {}
for line in open("projects_2016_04_28_08_30_42.tsv"):
    row = line.rstrip("\n").split('\t')
    proj = row[0]
    tumorType = row[1].split("-")[0].strip()
    projectInfo[proj] = {"tumor_type" :tumorType, "project_name" : row[1], "primary_site" : row[2]}

print "reading icd10 codes"
icd10ToText = {}
for line in open("/hive/data/outside/pubs/icd10/icd10cm_order_2016.txt"):
    #fs = line.rstrip("\n").split("\t")
    code = line[6:14].strip()
    desc = line[16:77].strip()
    icd10ToText[code] = desc

mutFnames = glob.glob("byProject/simple_somatic_mutation.open.*.tsv.gz")
#mutFnames = glob.glob("byProject/simple_somatic_mutation.open.ALL-US.tsv.gz")

vcfFh = None
projMetaDict = {}

lastProjCode = None
lastSampleId = None

for mutFname in mutFnames:
    print "reading donor info"
    donorInfo = {}
    donorFname = mutFname.replace("simple_somatic_mutation.open", "donor")
    for row in nextRow(gzip.open(donorFname)):
        #fs = line.rstrip("\n").split("\t")
        donorDict = row._asdict()
        diagText = icd10ToText.get(row.donor_diagnosis_icd10.replace(".", ""), "")
        donorDict["diagnosis"] = diagText
        donorInfo[row.icgc_donor_id] = donorDict

    print "reading specimen info"
    specFname = mutFname.replace("simple_somatic_mutation.open", "specimen")
    specInfo = {}
    for row in nextRow(gzip.open(specFname)):
        #fs = line.rstrip("\n").split("\t")
        specInfo[row.icgc_specimen_id] = row

    proj = basename(mutFname).split('.')[2]

    print "Reading", mutFname
    for row in nextRow(gzip.open(mutFname)):
        # if needed, flush the meta data
        if row.project_code!=lastProjCode and lastProjCode is not None:
            writeMeta(lastProjCode, projMetaDict)
            projMetaDict = {}

        # if needed, open a new vcf file handle
        if row.icgc_sample_id!=lastSampleId:
            #vcfFh = openVcf(row)
            if vcfFh is not None:
                vcfFh.close()
                cmd = ["bgzip", "-f", vcfFh.name]
                assert(subprocess.call(cmd)==0)

            # create a new sample meta info dict for every new sample id
            sampleData = {}
            # take donor data 
            donorDict = donorInfo[row.icgc_donor_id]
            sampleData.update(donorDict)
            # merge in specimen data
            specDict = specInfo[row.icgc_specimen_id]._asdict()
            sampleData.update(specDict)
            # merge in project data
            sampleData.update(projectInfo[row.project_code])

            vcfFname = row.project_code+"__"+row.submitted_sample_id+".vcf"
            sampleData["vcf_filename"] = vcfFname+".gz"
            sampleData["reference"] = "GRCh37"
            sampleData["sample_id"] = row.submitted_sample_id
            sampleData["icgc_sample_id"] = row.icgc_sample_id

            # pull experiment pipeline/platform info from mutation into sample
            rowDict = row._asdict()
            assert(row._fields[30:]==expFields)
            for field in expFields:
                sampleData[field] = rowDict[field]

            # remove keys with empty vals
            sampleDataNew = {}
            for key, val in sampleData.iteritems():
                if val!="":
                    sampleDataNew[key] = val
            sampleData = sampleDataNew

            projMetaDict[row.submitted_sample_id] = sampleData

            vcfPath = "ledgerSubmit/"+vcfFname
            print "opening new vcf file %s" % vcfPath
            vcfFh = open(vcfPath, "w")
            vcfFh.write(vcfHead)

        # a few basic checks
        int(row.chromosome_start)
        int(row.chromosome_end)
        assert(row.assembly_version=="GRCh37")
        if not (row.chromosome in validChroms):
            print row.chromosome
            assert(False)

        infoParts = []
        if row.mutated_from_allele != row.reference_genome_allele:
           infoParts = ["mutated_from_allele=%s" % row.mutated_from_allele ]
        if row.total_read_count!="":
            infoParts.append("TC=%s" % row.total_read_count)
        if row.mutant_allele_read_count!="":
            infoParts.append("MC=%s" % row.mutant_allele_read_count)
        if row.gene_affected!="":
            infoParts.append("GA=%s" % row.gene_affected)
        if row.transcript_affected!="":
            infoParts.append("TA=%s" % row.transcript_affected)
        if row.verification_status!="":
            infoParts.append("VS=%s" % row.verification_status)
        if row.verification_platform!="":
            infoParts.append("VP=%s" % row.verification_platform)
        infoField = ";".join(infoParts)

        #eight fields: CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
        vcfRow = [row.chromosome, row.chromosome_start, row.icgc_mutation_id, row.reference_genome_allele, row.mutated_to_allele, ".", ".", infoField]
        vcfFh.write("\t".join(vcfRow))
        vcfFh.write("\n")

        lastSampleId = row.icgc_sample_id
        lastProjCode = row.project_code

writeMeta(lastProjCode, projMetaDict)
