# functions to deal with genomics datasets, bed files etc.
# mostly wrappers around kent src tools

import collections, logging, maxTables, types
from os.path import *

# GLOBAL CONFIG VARS
genomeDataDir = "/hive/data/genomes"

# STRUCTS
Bed3 = collections.namedtuple("bed3", ["chrom", "start", "end"])

def listToStr(list, sep="\t"):
    " cast list to strings and join into tab sep string "
    list = [str(x) for x in list]
    return "\t".join(list)

def listToStrLn(list, sep="\t"):
    return toString(list)+"\n"

def getChromSizes(db):
    " return chrom sizes as dict chromName -> size "
    if not isfile(db):
        chromSizesFname = join(genomeDataDir, db, "chrom.sizes")
    else:
        chromSizesFname = db
    chromSizes = maxTables.parseDict(chromSizesFname, valueType=types.IntType)
    return chromSizes

def chromSplit(db, step):
    """ 
    yield bed3 features with chromosomes split into 
    nonoverlapping pieces of size size 
    """
    chromSizes = getChromSizes(db)
    for chrom, size in chromSizes.iteritems():
        size = int(size)
        for start in range(0, size, step):
            end = start + step -1
            if end > size:
                end = size
            yield Bed3(chrom, start, end)

