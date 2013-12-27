from optparse import OptionParser
import cgitb
import cgi
import sys
import os.path
import tabfile

# generic python library for processing input from the command line OR html webpages
# - parameters are specified by their properties (file, string, int, etc)
# - html web forms or a command line interface is generated from these
# - goal: add the extension .cgi to your scripts and they become simple webpages WITHOUT changing anything else!
# - why: urls can be sent to your biologist "friends" :-) but scripts cannot


# Example:
# ip = InputProcessor("ImageSorter", "Generate sorted reference sheets from Ghost/Aniseed in-situ pictures")
#p = Parameter("textLinesArg", "genes", desc="list of genes, optional: annotation, separated by tab; optional: in-situ id, separated by another tab")
#ip.add_parameter(p)
#p = Parameter("textLines", "aniseedAnnot", shortOption="-a", desc="default expression annotation", cgiHide=True, cgiDefault="aniseed-etb-mtb.lst", textFormat="dictList")
#ip.add_parameter(p)
# input = ip.get_input()

class Parameter:
    """ a parameter to a script """
    def __init__(self, ptype, name, textFormat="", shortOption="", longOption="", desc="", longDesc="", cgiDefault="", clDefault="", cgiHide=False, required=False):
        """ ptype is one of 
             - "string" (a single value)
             - "boolean" (yes/no) 
             - "textLines" (textarea/option for textfile), 
             - "textLinesArg" (textarea/argument), 
            textFormat is one of dictList, dict """
        self.ptype=ptype.lower()
        if ptype=="textlines" or ptype=="textlinesarg":
            self.clType="string"
        else:
            self.clType="string"
        self.name=name
        self.desc=desc
        self.longDesc=longDesc
        self.cgiDefault=cgiDefault
        self.clDefault=clDefault
        self.cgiHide=cgiHide
        self.shortOption=shortOption
        self.longOption=longOption
        self.textFormat=textFormat.lower()
        self.required=required

    def __repr__(self):
        lines = []
        lines.append("Parameter-Object")
        lines.append("ptype="+self.ptype)
        lines.append("name="+self.name)
        lines.append("desc="+self.desc)
        lines.append("longDesc="+self.longDesc)
        lines.append("cgiDefault="+self.cgiDefault)
        lines.append("clDefault="+self.clDefault)
        lines.append("cgiHide="+str(self.cgiHide))
        return ", ".join(lines)


class InputProcessor:
    def __init__(self, title, clDesc="", cgiDesc="", inputFormat="", cgiMethod="post"):
        """ scriptName is the name of .cgi-file """
        self.params = []
        self.scriptName = os.path.split(sys.argv[0])[1]
        self.cgiMode = False
        self.cgiMethod = cgiMethod
        self.title = title
        self.clDesc = clDesc
        self.cgiDesc = cgiDesc
        self.inputFormat = inputFormat
        if cgiDesc=="":
            self.cgiDesc = clDesc
        elif clDesc=="":
            self.clDesc = cgiDesc

    def add_parameter(self, p):
        self.params.append(p)

    def print_form(self):
        """ print html form to stdout """
        print '<form action="%s" method="%s">' % (self.scriptName, self.cgiMethod)
        for p in self.params:
            if p.cgiHide:
                continue
            elif p.ptype=="textlines" or p.ptype=="textlinesarg":
                print '<textarea name="%s" cols="60" rows="20">%s</textarea>' % (p.name, p.cgiDefault)
            elif p.ptype=="boolean":
                checked = ""
                if p.cgiDefault.lower()=="true":
                    checked="CHECKED"
                print '<input type="CHECKBOX" name="%s" %s>%s</input>' % (p.name, checked, p.desc)
            elif p.ptype=="string":
                print '%s: <input type="text" size="20" name="%s" value="%s"></input>' % (p.desc, p.name, p.cgiDefault)

            else:
                print "internal error: illegal ptype=%s" % p.ptype
            
            if p.longDesc:
                print '<br><small>%s</small>' % p.longDesc

            print '<p>'

        print '<input TYPE="submit" NAME="submit" >     '
        #print '<input TYPE="submit" NAME="submit" VALUE="     OK     ">     '
        print '<input TYPE="reset" NAME="reset" VALUE="   Reset "><br>'
        print '</form>'

    def print_css(self):
        print """
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
        .geneHead { width=200px; font-size:7px; font-weight:bold; background-color: #DDDDDD; border-top: 1px solid #BBB; margin-top: 2em; }
        .geneDetails { width=200px; border-bottom: 1px solid #BBB; }
        .geneGood    { background-color: green; }
        .geneNotGood { background-color: red; }
        .geneDetails { width=200px; border-bottom: 1px solid #BBB; font-weight: bold; }
        .regions     { font-size: 7px; }
        .regionHead  { font-size: 7px; font-weight:bold; background-color: #EEEEEE; }
        .exprText    { font-size:7px; font-weight: normal; text-align:left; padding-left: 10px; vertical-align: top; }
        </style>
        """

    def print_html(self):
        """ print complete html webpage with form """
        print "<html>"
        print "<head>"
        print "<title>%s</title" % self.title
        print "</head>"
        print "<body>"
        print "<h2>%s</h2>" % (self.title)
        print "<h4>%s</h4>" % (self.cgiDesc+":")
        print "%s<br>" % (self.inputFormat)
        self.print_form()
        print "</body>"
        print "</html>"

    def generate_parser(self):
        arguments = [p.name for p in self.params if (p.ptype=="textlinesarg" or p.required==True)]
        parser = OptionParser("usage: %s %s [options] - %s" % (self.scriptName, " ".join(arguments), self.clDesc) )
        for p in self.params:
            if p.ptype=="textlinesarg" or p.required:
                continue
            elif p.ptype=="boolean":
                parser.add_option(p.shortOption, p.longOption, dest=p.name, action="store_true", help=p.desc, default=p.clDefault)
            else:
                parser.add_option(p.shortOption, p.longOption, dest=p.name, action="store", help=p.desc, type=p.clType, default=p.clDefault)
        return parser

    def cgi_params_ok(self):
        """ check if all parameters are present """
        form = cgi.FieldStorage()
        error = False
        for p in self.params:
            if p.cgiHide or not p.required:
                continue
            if form.getfirst(p.name)==None and (p.required and not p.ptype=="boolean"):
                return form, False
        #if form.getfirst("submit"): # modified to allow non-submit requests
        if len(form)!=0:
            return form, True
        else:
            return form, False

    def cgi_results(self, form):
        values = {}
        for p in self.params:
            if (not form.has_key(p.name)) and p.cgiHide:
                val = p.cgiDefault
            else:
                val = form.getfirst(p.name)

            # postproc
            if p.ptype=="boolean":
                if val==None:
                    val = False
                else:
                    val = True
            elif p.ptype=="textlines" or p.ptype=="textlinesarg":
                val = self.paramFile(p, val, True)

            # store
            values[p.name]=val
        return values

    def paramFile(self, param, value, cgiMode=False):
        """ return value OR contents of value (name-of-file) if conditions are met to read in the file """
        if value!="":
            if (param.ptype.lower()=="textlinesarg" and not cgiMode) or param.ptype.lower()=="textlines":
                if param.textFormat=="":
                    return tabfile.slurplist(value)
                elif param.textFormat=="dictlist":
                    return tabfile.slurpdictlist(value)
                elif param.textFormat=="dict":
                    return tabfile.slurpdict(value)
            else:
                return value
        else:
            return value

    
    def cl_results(self):
        parser = self.generate_parser()
        options, args = parser.parse_args()
        values = {}
        optDict = options.__dict__

        # read options, open and parse files
        i=0
        for p in self.params:
            if p.ptype=="textlinesarg":
                values[p.name] = self.paramFile(p, args[i])
            elif p.required and p.ptype=="string":
                values[p.name] = args[i]
            else:
                values[p.name] = self.paramFile(p, optDict[p.name])
            i+=1

        return values

    def get_input(self, charset="utf-8"):
        if self.scriptName.endswith(".cgi"):
            cgitb.enable(); # replace python handler with html one
            form, cgiOK = self.cgi_params_ok()
            print "Content-Type: text/html; charset=%s" % charset    
            print                               # blank line, end of headers
            if cgiOK:
                self.cgiMode=True
                return self.cgi_results(form)
            else:
                self.print_html()
                sys.exit(0)
        else:
            if len(sys.argv)==1:
                parser = self.generate_parser()
                parser.print_help()
                sys.exit(1)
            else:
                self.cgiMode=False
                return self.cl_results()



    



