# these are config options that mostly relate to 
# text2genome's webinterface

# search for python modules in the 'lib/'-directory

import sys, os, ConfigParser, logging

config = ConfigParser.ConfigParser()
config.read(['t2g.conf', os.path.expanduser('~/.t2g.conf')])

section = None

def initFromString(section, strings):
    """ converts string like 'bla=test test=2 bla=4' to a dict and initilizes
    config from it. Used to update config from command line parameters """
    for s in strings:
        key, value = s.split("=")
        config.set(section, key, value)

def setSection(sec):
    global section
    section = sec

def getValue(key, default):
    return get(section, key, default)

def getPath(key, default=None):
    value = getValue(key, default)
    if value!=None:
        return os.path.expanduser(value)
    else:
        return None

def get(section, key, default):
    """ retrieve value for key in config file, if not found return default """

    if section in config.sections():
        try:
            val = config.get(section, key)
            logging.info("Config %s/%s=%s" %(section, key, val))
            return val
        except ConfigParser.NoOptionError:
            logging.info("Config value %s/%s not found: returning %s" %(section, key, default))
            return default
    else:
        logging.info("Config %s/%s: section not found, returning %s" %(section, key, default))
        return default

def getBool(section, key, default):
    return bool(int(get(section, key, default)))

def getAllPrefix(section, prefix):
    """ get all values in section that start with prefix """
    prefix=prefix.lower()
    values = []
    for name, value in config.items(section):
        if name.startswith(prefix):
            key = name.replace(prefix+".", "")
            values.append((key, value))
    logging.info("Config %s/prefix %s: %s" %(section, prefix, str(values)))
    return values

def sqlConnStringToDict(connString):
    """ convert a string in format host:port,username,passwd,db to a dictionary that can be passed directly to mysqldb.connect """
    if connString==None:
        return {}

    host,user,passwd,db = connString.split(",")

    port = None
    fs = host.split(":")
    host = fs[0]
    if len(fs)>1:
        port=int(fs[1])

    dict = {"host":host, "user":user, "passwd":passwd, "db":db, "read_default_file":"/etc/my.cnf"}
    if port:
        dict["port"]=port
    
    return dict

#
mysqlHost = "127.0.0.1"
mysqlPort = 3306
mysqlUser = "mhaeussl"
mysqlPassword = ""
mysqlDb = "text2genome"
menu = [ 
                ("Text2Genome" , [ 
                    ("About","about.cgi"), 
                    ("Search","search.cgi"), 
                    ("Browse","browserOverlay.cgi"),
                    ("Download","download.cgi"),
                    ("API","api.cgi")
                    ] ) ,
                ("About us" , [ ("Bergman Lab", "http://www.bioinf.manchester.ac.uk/bergman/") ] )
            ]

host = "http://max.smith.man.ac.uk"
baseDir = host+"/t2g"
# We need the FULL URL of the inspector.cgi script on your webserver
detailsUrl= baseDir+'/inspector.cgi'
# We need the FULL URL of the DAS root on your webserver
dasUrl = baseDir+'/das'
# We need the FULL URL of the ucsc tracks on your webserver
bedDir= baseDir+'/bed/'

ensemblServer = "sep2009.archive.ensembl.org"

def dasUri(genomeId):
    """ used by the das server and the viewer to retrieve the DAS uri """
    return "t2g_ensembl56_taxId%d" % (genomeId)
