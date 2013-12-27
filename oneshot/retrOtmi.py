#!/usr/bin/env python
import util
import simplexmlparse
import os

""" script to download otmi-files from nature into subdirectory files/ """
""" uses simplexmlparse """

TEMPLATE= """<?xml version="1.0" encoding="UTF-8"?>
<opml version="1.0" xmlns:simplexmlparse="http://evanjones.ca/simplexmlparse">
  <head>
    <title>OTMI: Journals</title>
    <dateCreated>Fri, 15 Dec 2006 11:00:00 +0000</dateCreated>
    <dateModified>Thu, 1 Mar 2007 11:37:00 +0000</dateModified>
    <ownerName>OTMI (Open Text Mining Interface)</ownerName>
    <ownerEmail>otmi@nature.com</ownerEmail>
  </head>
  <body>
    <outline simplexmlparse:count="*" type="opml" text="The Pharmacogenomics Journal" opmlUri="http://www.nature.com/tpj/otmi/tpj.opml"/>
  </body>
</opml>
"""

template2 = """<?xml version="1.0" encoding="UTF-8"?>
<opml version="1.0" xmlns:simplexmlparse="http://evanjones.ca/simplexmlparse">
  <head>
    <title>OTMI: Nature</title>
    <dateCreated>Fri, 15 Dec 2006 11:00:00 +0000</dateCreated>
    <dateModified>Fri, 15 Dec 2006 11:00:00 +0000</dateModified>
    <ownerName>OTMI (Open Text Mining Interface)</ownerName>
    <ownerEmail>otmi@nature.com</ownerEmail>
  </head>
  <body>
      <outerOutline text="test">
          <outline simplexmlparse:count="*"  type="opml" text="Nature 433(7021)" opmlUri="http://www.nature.com/nature/journal/v433/n7021/otmi/otmi-manifest.opml" gzipUri="http://www.nature.com/nature/journal/v433/n7021/otmi/otmi-contents.tar.gz" />
      </outerOutline>
  </body>
</opml>"""

# ---- MAIN ------
url = "http://www.nature.com/otmi/journals.opml"
opml = util.parseXmlUrl(TEMPLATE, url)
for outline in opml.body.outline:
    url= outline.opmlUri
    opmlStr = util.httpGet(url).read()
    # hack! workaround for simplexmlparse's limitations
    # this can easily break but will trigger an exception if it does one day
    opmlStr = opmlStr.replace("<outline text=", "<outerOutline text=").replace("</outline>", "</outerOutline>")
    opmlStr = opmlStr.replace("<outline>", "<outerOutline>").replace("</outline>", "</outerOutline>")
    journalOpml = util.parseXml(template2, opmlStr)
    for o2 in journalOpml.body.outerOutline.outline:
        gzipUrl = o2.gzipUri
        print "downloading "+gzipUrl
        util.httpDownload(gzipUrl, "files/file.tar.gz")
        os.system("cd files; tar xvfz file.tar.gz")
