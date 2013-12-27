#!/usr/bin/env python
import logging, re, urllib, StringIO
try:
    import xml.etree.cElementTree as etree
except:
    try:
        import elementtree.ElementTree as etree
    except:
        logging.error("Cannot load any elementtree library, sure you're running python>=2.5?")
    

class XmlParser:
    """ class to represent an xml tree (using ElementTree)
        Functions Accept PATH which is a /a/b/c style xpath-like expression to refer to elements
        PATH is not a complete XPATH implementation

        getText... functions return just a string
        getXml... functions return an XmlParser-object
        ...First  functions get only the first instance
        ...All    functions return an iterator
    
    >>> xp = XmlParser(string="<fruit><apple size='big'>boskoop</apple><apple size='small'>granny smith</apple><pear>mypear</pear></fruit>") 
    >>> xp.getTextFirst("pineapple", default="NothingAtAll")
    'NothingAtAll'
    >>> xp.getTextFirst("apple")
    'boskoop'
    >>> list(xp.getTextAll("apple"))
    ['boskoop', 'granny smith']
    >>> list(xp.getTextAll("apple", reqAttrDict={'size':'big'}))
    ['boskoop']
    
    """
    def __init__(self, string=None, url=None, root=None):
        self.root=None
        if string:
            self.fromString(string)
        elif url:
            self.fromUrl(url)
        elif root:
            self.root=root

    def getText(self):
        return self.root.text

    def fromString(self, string, removeNamespaces=False):
        root = etree.fromstring(string)
        if removeNamespaces:
            logging.debug("Stripping all namespaces")
            stripNamespaces(root)
        self.root = root

    def fromUrl(self, url, removeNamespaces=False, stopWords=[]):
        logging.debug("Retrieving %s" % url)
        text = urllib.urlopen(url).read()
        self.fromString(text, removeNamespaces=removeNamespaces)
        #for w in stopWords:
            #if w in text:
                #return None

    def _removeNamespaces(self):
        """ removes all namespaces from elementtree IN PLACE """
        root = self.root
        for el in root.getiterator():
            if el.tag[0] == '{':
                el.tag = el.tag.split('}', 1)[1]

    def _hasAttribs(self, el, reqAttrDict):
        for attr, value in reqAttrDict.iteritems():
            if el.attrib.get(attr, None)!=value:
                return False
        return True

    def getTextFirst(self, path, reqAttrDict=None, default=None):
        """ return text between elements given path 
            reqAttrDict is in the format attrName -> value
        """
        xml = self._getElFirst(path, reqAttrDict)
        if xml != None:
            return xml.text
        else:
            return default
        
    def getTextAll(self, path, reqAttrDict=None):
        for el in self._getElAll(path, reqAttrDict):
            yield el.text

    def _getElFirst(self, path, reqAttrDict):
        found = False
        for el in self._getElAll(path, reqAttrDict):
            found = True
            return el
        return None

    def _getElAll(self, path, reqAttrDict):
        found = False
        elIter = self.root.findall(path)
        for el in elIter:
            if reqAttrDict == None or self._hasAttribs(el, reqAttrDict):
                found = True
                yield el
        
    def getXmlFirst(self, path, reqAttrDict=None, default=None):
        el = self._getElFirst(path, reqAttrDict)
        if el==None:
            return default
        else:
            return XmlParser(root=el)

    def getXmlAll(self, path, reqAttrDict=None):
        for el in self._getElAll(path, reqAttrDict):
            yield XmlParser(root=el)

# =============== XML WRITER ========================
def makeAttrStr(dict):
    """ given a dict {"font":"fat", "my_Type":"big"} create string 'font="fat" my-type="big"' """
    if len(dict)==0:
        return ""

    attrList=[]
    for key, val in dict.iteritems():
        partString = key.replace("_", "-") +'="%s"' % str(val)
        attrList.append(partString)
    return " "+" ".join(attrList)
            
class XmlWriter:
    """ a class that makes writing xmls files easier 
    example:
    x = xmlWriter(sys.stdout)
    x.tagStart("a", href="http://www.sgi.com")
    x.tagEnd("a")
    x.tag("a", content="link", href="http://www.altavista.com")
    """

    def __init__(self, out):
        if isinstance(out, file):
            self.f = out
        else:
            self.f = open(out, "w")

    def writeLn(self, str):
        self.f.write(str)
        self.f.write("\n")

    def tag(self, tagName, content=None, **attr):
        if content==None:
            self.writeLn('<%s%s />' % (tagName, makeAttrStr(attr)))
        else:
            self.writeLn('<%s%s>%s</%s>' % (tagName, makeAttrStr(attr), content, tagName))

    def tagStart(self, tagName, **attr):
        self.writeLn('<%s%s>' % (tagName, makeAttrStr(attr)))

    def tagEnd(self, tagName, **attr):
        self.writeLn('</%s%s>' % (tagName, makeAttrStr(attr)))

    ## ====== SVG TAGS ===================
    def svgStart(self, **attr):
        # <svg version="1.1" baseProfile="full" > 
        attr["baseProfile"]="full"
        attr["xmlns"]="http://www.w3.org/2000/svg"
        attr["version"]="1.1"
        self.tagStart("svg", **attr)

    def svgEnd(self, **attr):
        self.tagEnd("svg")

    def text(self, text, **attr):
        self.tag("text", text, **attr)

    def gStart(self, **attr):
        self.tagStart("g", **attr)

    def gEnd(self, **attr):
        self.tagEnd("g")

    def rect(self, **attr):
        self.tag("rect", **attr)

    def line(self, **attr):
        self.tag("line", **attr)

    def transform(self, **attr):
        self.tag("transform", **attr)

    def stroke(self, **attr):
        self.tag("stroke", **attr)

    def path(self, **attr):
        self.tag("path", **attr)

# === BUFFER BUSINESS ===== 
# writing can be temporarily redirected to a buffer
# this buffer can then later be written
    def writeToBuffer(self):
        " will direct all subsequent writing to temporary string instead of file "
        self.buffer = StringIO.StringIO()
        self.realFile = self.f
        self.f = self.buffer

    def writeBuffer(self):
        " write out buffer string"
        self.f.write(self.buffer.getvalue())

    def writeToFile(self):
        " switch back to direct file writing "
        self.f=self.realFile

    def clearBuffer(self):
        self.buffer = StringIO.StringIO()

# ----- 
if __name__ == "__main__":
    import doctest
    doctest.testmod()
