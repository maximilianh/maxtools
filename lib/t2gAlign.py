from sys import *
import logging, os
import util
logger = logging.getLogger("t2gAlign")

# copied from http://code.activestate.com/recipes/302594/ 
def relpath(target, base=os.curdir):
    """
    Return a relative path to the target from either the current dir or an optional base dir.
    Base can be a directory specified either as absolute or relative to current dir.
    """

    if not os.path.exists(target):
        raise OSError, 'Target does not exist: '+target

    if not os.path.isdir(base):
        raise OSError, 'Base is not a directory or does not exist: '+base

    base_list = (os.path.abspath(base)).split(os.sep)
    target_list = (os.path.abspath(target)).split(os.sep)

    # On the windows platform the target may be on a completely different drive from the base.
    if os.name in ['nt','dos','os2'] and base_list[0] <> target_list[0]:
        raise OSError, 'Target is on a different drive to base. Target: '+target_list[0].upper()+', base: '+base_list[0].upper()

    # Starting from the filepath root, work out how much of the filepath is
    # shared by base and target.
    for i in range(min(len(base_list), len(target_list))):
        if base_list[i] <> target_list[i]: break
    else:
        # If we broke out of the loop, i is pointing to the first differing path elements.
        # If we didn't break out of the loop, i is pointing to identical path elements.
        # Increment i so that in all cases it points to the first differing path elements.
        i+=1

    rel_list = [os.pardir] * (len(base_list)-i) + target_list[i:]
    return os.path.join(*rel_list)

def getTargetFiles(dir, ext):
    fileData = os.walk(dir)
    files = []
    for dirpath, dirnames, filenames in fileData:
        for file in filenames:
            if file.endswith(ext):
                files.append(os.path.join(dirpath, file))
    return files

def filterTargetFiles(dbFiles, queryFile):
    """ keep only dbFiles-entries that contain queryFile somewhere in their paths """
    queryGenome=queryFile.replace(".fa", "")
    queryGenome=queryGenome.lower()

    newDbFiles = []
    for dbPath in dbFiles:
        dbDirs = dbPath.split("/")
        dbDirs = [db.lower() for db in dbDirs]
        if queryGenome in dbDirs:
            newDbFiles.append(dbPath)

    return newDbFiles

def runAligner(blastPath, blastType, inputFilename, targetPath,
    targetExtension, targetFilterWord, outDir, test, aligner,
    flatOutDir, clusterSystem, trySub, queue, logPath, addOptions):

    bsub=False
    qsub=False
    if clusterSystem!=None:
        if clusterSystem.lower()=="sge":
            qsub=True
        elif clusterSystem.lower()=="lsf":
            bsub=True

    if not os.path.isdir(outDir):
        logger.error("output directory %s does not exist" % (outDir))
        exit(1)
    if not os.path.isdir(targetPath):
        logger.error("target database directory %s does not exist" % (targetPath))
        exit(1)

    dbFiles = getTargetFiles(targetPath, targetExtension)

    if targetFilterWord:
        dbFiles = filterTargetFiles(dbFiles, targetFilterWord)

    if len(dbFiles)==0:
        logger.error("error: No blast/blat files (*.nin or *.2bit) found for input file %s in subdirectories of %s!\n" % (inputFilename, targetPath))
        exit(1)
    logger.info("Found %d files to BLAST/BLAT on\n" % len(dbFiles))

    subCount=0

    if aligner=="bwa":
        logging.info("Splitting fasta files in long (>200bp) and short (<200bp)")

        shortFastaFilename = inputFilename+".short"
        longFastaFilename  = inputFilename+".long"

        util.execCmdLine("faFilter -maxSize=200 %s %s" % (inputFilename, shortFastaFilename))
        util.execCmdLine("faFilter -minSize=200 %s %s" % (inputFilename, longFastaFilename))

    for db in dbFiles:

        if flatOutDir:
            # place all output into one directory, requires that filename includes genome
            outfile = os.path.join(outDir, os.path.basename(db))
        else:
            # write to same relative directory from outDir, create all necessary directories
            outfile = os.path.join(outDir, relpath(db, base=targetPath))
            if not test:
                dirname = os.path.dirname(outfile)
                if not os.path.isdir(dirname):
                    os.makedirs(dirname)

        outfile = outfile.replace(targetExtension,"").replace(".fa","")

        infile=db
        infile=os.path.normpath(infile)
        optionString = addOptions

        # CONSTRUCT COMMAND LINES
        # BLAST 
        if aligner=="blast":
            db = db.replace(targetExtension,"")
            optionString += " -m 8"
            if blastPath==None or len(blastPath)==0:
                blastExec="blastall"
            else:
                blastExec=blastPath
            cmdLine = "%s -p %s %s -i %s -d %s" % (blastExec, blastType, optionString, inputFilename, db)
            cmdLine += " -o %s" % (outfile+".blast")

        # BLAT 
        elif aligner=="blat":
            cmdLine = "blat -noHead %s %s %s %s" % (infile, inputFilename, (outfile+".psl"), optionString)

        # BWA 
        elif aligner=="bwa":
            db = db.replace(targetExtension,"")
            dbStem = os.path.basename(db)
            shortSaiFilename = outfile+".sai"
            shortSamFilename = outfile+".short.sam"
            longSamFilename = outfile+".long.sam"
            outfile = outfile+".sam"
            
            cmdLine = "bwa aln %s %s > %s; " % (db, shortFastaFilename, shortSaiFilename)
            cmdLine +=  "bwa samse %s %s %s > %s; " % (db, shortSaiFilename, shortFastaFilename, shortSamFilename)
            cmdLine +=  "bwa bwasw %s %s > %s; " % (db, longFastaFilename, longSamFilename)
            cmdLine +=  "cat %s %s > %s; " % (shortSamFilename, longSamFilename, outfile)
            cmdLine +=  "rm -f %s %s %s" % (shortSamFilename, longSamFilename, shortSaiFilename)
        else:
            assert(False) # illegal aligner?

        # LSF cluster scheduler
        if bsub:
            cmdLine = "bsub -q "+queue+" '"+cmdLine+"'"
        # SUN GRID ENGINE cluster scheduler
        elif qsub:
            logFileStr = ""
            if logPath:
                logFullPath = os.path.join(logPath,"t2g-$JOB_ID-$JOB_NAME.out")
                logFileStr = "-o '"+logFullPath+ "' "

            queueStr = "-q "+queue # will prevent stopping of jobs > 24hrs

            jobName = os.path.splitext(os.path.basename(infile))[0]
            cmdLine = 'qsub %s -N %s -cwd -V -j y %s -b y "' % (logFileStr, jobName, queueStr) +cmdLine+ '"' 

        logger.info( "Running: "+cmdLine)

        if not test:
            ret = os.system(cmdLine)
        else:
            print cmdLine
            ret = 0
        subCount+=1

        if trySub and (subCount>=3):
            logger.info( "--try was set, stopped after %d submissions" % subCount)
            break

        if ret != 0:
            logger.info("error occured, return code != 0\n")
            exit(1)
        
    if test:
        logger.info("Test mode activated, no command was run, remove testMode " +
        "to actually run something.")
    else:
        logging.info("BLAST jobs completed")

if __name__ == "__main__":
    import doctest
    doctest.testmod()

