#!/usr/bin/python
#
"""A simple framework for validating XML documents and parsing them into
Python objects. Consider it like a VERY simple version of RelaxNG. More
documentation can be found at http://evanjones.ca/software/simplexmlparse.html

This module has a main routine which prints the Python object tree created by
parsing the document with the template."""

__author__ = 'Evan Jones <ejones@uwaterloo.ca>'
__url__ = 'http://evanjones.ca/software/simplexmlparse.html'
__date__ = '$Date: 2006/04/16 16:58:06 $'.split()[1].replace('/', '-')
__version__ = '$Revision: 1.17 $'

import xml.sax as sax
from cStringIO import StringIO

SIMPLEPARSE_NAMESPACE = "http://evanjones.ca/simplexmlparse"
SIMPLEPARSE_REFERENCE = (SIMPLEPARSE_NAMESPACE, 'ref')
SIMPLEPARSE_COUNT = (SIMPLEPARSE_NAMESPACE, 'count')

# Defines the attributes that are permitted in SIMPLEPARSE_NAMESPACE
SIMPLEPARSE_ATTRIBUTES = {
	SIMPLEPARSE_REFERENCE: 1,
	SIMPLEPARSE_COUNT: 1,
}


E_INVALID_CHILD = 2
E_INVALID_DOCTAG = 3
E_NO_CHILDREN = 5
E_MISSING = 7
E_TEXT = 8
E_MULTIPLE = 9
E_COUNT = 10
E_INVALID_ATTRIBUTE = 11
E_CONFLICT = 12
E_MISSING_ATTRIBUTE = 13
E_INVALID_IDENTIFIER = 14
# Template parse errors
E_ELEMENT_CONFLICT = 15
E_UNDEFINED_REF = 16
E_REF_NOT_EMPTY = 17
E_REF_HAS_ATTRIBUTE = 18

ERROR_CODES = {
	E_INVALID_CHILD: "Element '%s' is not permitted to have child element '%s'",
	E_INVALID_DOCTAG: "Invalid document tag '%s' (must be '%s')",
	E_NO_CHILDREN: "Element '%s' cannot contain any child elements (child tag: '%s')",
	E_MISSING: "Element '%s' is missing a required element: '%s'",
	E_TEXT: "Element '%s' cannot contain non-whitespace text%s",
	E_MULTIPLE: "Element '%s' cannot contain multiple '%s' elements",
	E_COUNT: "Element '%s' has an invalid simplexmlparse:count attribute value: '%s' (must be '*', '?', '1', or '+')",
	E_INVALID_ATTRIBUTE: "Element '%s' is not permitted to have attribute: '%s'",
	E_CONFLICT: "Element '%s' has conflicting attribute and child both named '%s'",
	E_MISSING_ATTRIBUTE: "Element '%s' is missing a required attribute: '%s'",
	E_INVALID_IDENTIFIER: "Invalid identifier: '%s' (%s)",
	E_ELEMENT_CONFLICT: "Element name '%s' has already been defined (elements must be unique: use a reference)%s",
	E_UNDEFINED_REF: "Reference from '%s' to undefined element '%s' (define it before the reference)",
	E_REF_NOT_EMPTY: "Reference from '%s' to '%s' not empty (references MUST be empty)",
	E_REF_HAS_ATTRIBUTE: "Reference '%s' has attribute '%s' (references MUST not have attributes: put them on the original element)",
	}

def validateIdentifier( namespaceName, locator ):
	'''Internal function. Returns the qualified name converted to a valid Python
	identifier. This will raise the appropriate ParseError if the identifier is not valid.'''
	
	try:
		identifier = namespaceName[1].encode( 'ascii' )
	except UnicodeEncodeError:
		raise ParseError( E_INVALID_IDENTIFIER, locator, namespaceName, 'not ASCII' )
	
	if identifier.startswith( '_' ):
		raise ParseError( E_INVALID_IDENTIFIER, locator, namespaceName, 'leading underscores are not permitted' )
	
	if ':' in identifier:
		raise ParseError( E_INVALID_IDENTIFIER, locator, namespaceName, 'colons are not permitted' )
	
	if '-' in identifier:
		raise ParseError( E_INVALID_IDENTIFIER, locator, namespaceName, 'dashes are not permitted' )
	
	if '.' in identifier:
		raise ParseError( E_INVALID_IDENTIFIER, locator, namespaceName, 'periods are not permitted' )
	
	return identifier

class ParseError( Exception ):
	'''Exception raised by the SimpleXMLParse module. This type has the following
	member variables:
	
	code: A numeric code representing the reason for the error. All the codes are
		defined in this module, prefixed by "E_"
	message: A string representing the reason for the error.
	line: The line number where the error occurred.
	column: The column number where the error occurred.'''
	
	__slots__ = ( 'code', 'message', 'line', 'column', )
	
	def __init__( self, code, locator, tagOne, tagTwo ):
		'''Creates a new ParseError with code at XML location locator. The
		objects tagOne and tagTwo will be formatted into the message.'''
		
		assert( code in ERROR_CODES )
		
		self.code = code
		self.line = locator.getLineNumber()
		self.column = locator.getColumnNumber()
		
		t1 = tagOne[1]
		if tagOne[0] != None:
			t1 = tagOne[0] + ':' + t1
		
		if tagTwo == None:
			t2 = ''
		elif isinstance( tagTwo, basestring ):
			t2 = tagTwo
		else:
			t2 = tagTwo[1]
			if tagTwo[0] != None:
				t2 = tagTwo[0] + ':' + t2
		
		self.message = ERROR_CODES[code] % ( t1, t2 )
	
	def __str__( self ):
		'''Returns a string representation of this ParseError. This is
		equivalent to the "message" member variable, with the line and
		column numbers appended.'''
		
		m = self.message + " at line %d column %d" %( self.line, self.column )
		
		return m


class ElementBase( object ):
	'''The parser creates instances with this base class.'''
	
	__slots__ = ()
	
	def __init__( self, elementInfo ):
		'''Create a new instance using the element information in elementInfo.'''
		
		# If this instance should have text, initialize it to empty
		if elementInfo.containsText:
			self._text = ''
		
		# Make all the child nodes equal to None
		for name in elementInfo.attributes.iterkeys():
			# An attribute is None or a string
			setattr( self, name, None )
		
		# Make all the child nodes equal to None or empty lists
		for name in elementInfo.children.iterkeys():
			count = elementInfo.childrenCount[name]
			if count == COUNT_ZERO_OR_ONE or count == COUNT_ONE:
				# A single value is None or a string
				setattr( self, name, None )
			else:
				assert( count == COUNT_ZERO_OR_MORE or count == COUNT_ONE_OR_MORE )
				# Multiple values are a list
				setattr( self, name, [] )

COUNT_ZERO_OR_ONE = '?'
COUNT_ZERO_OR_MORE = '*'
COUNT_ONE = '1'
COUNT_ONE_OR_MORE = '+'

ELEMENT_COUNTS = ( COUNT_ZERO_OR_ONE, COUNT_ZERO_OR_MORE, COUNT_ONE, COUNT_ONE_OR_MORE )

class ElementInfo( object ):
	'''Stores the information about what is permitted by a specific
	element. This is built by the TemplateParser. After it has finished
	parsing the element, it calls createObjectType which creates the type
	for the object representing this element's data. Finally, the parser uses
	the ElementInfo instance to validate documents and create instances of
	the custom type.'''
	
	__slots__ = ( 'tagName', 'identifier', 'objtype', 'instance',
		'containsText', 'children', 'childrenCount', 'attributes' )
	
	def __init__( self, name, identifier ):
		'''Creates a new ElementInfo for the qualified element name, with
		the Python identifier.'''
		
		self.tagName = name
		self.identifier = identifier
		self.objtype = None
		self.instance = None
		
		# True if this element can contain non-whitespace text, False otherwise
		self.containsText = False
		
		# Maps identifiers to the ElementInfo for the child element
		self.children = {}
		# Maps identifiers to the counts for the child element
		self.childrenCount = {}
		# Maps identifiers to the count (either COUNT_ONE or COUNT_ZERO_OR_ONE) for the attribute
		self.attributes = {}
	
	def createObjectType( self ):
		'''After setting all the required attributes, this method will
		create a type to represent instances of this element.'''
		
		assert( self.objtype == None )
		
		# Permit attributes with the names of the permitted child elements
		slots = { '__slots__' : self.children.keys() }
		
		slots['__slots__'].extend( self.attributes.keys() )
		
		if self.containsText:
			slots['__slots__'].append( '_text' )
		
		self.objtype = type( self.identifier, (ElementBase,), slots )
	
	def create( self ):
		'''Creates a new instance of the object represented by this element.'''
		
		assert( self.instance == None )
		
		self.instance = self.objtype( self )
	
	def addChild( self, name, info, count ):
		'''Adds a permitted child element to this element.'''
		assert( self.objtype == None )
		assert( name not in self.children )
		assert( name not in self.attributes )
		assert( info != None )
		assert( count in ELEMENT_COUNTS )
		assert( not name.startswith( '_' ) )
		
		self.children[name] = info
		self.childrenCount[name] = count
	
	def addAttribute( self, identifier, count ):
		'''Adds a permitted attribute to this element.'''
		assert( self.objtype == None )
		assert( identifier not in self.children )
		assert( identifier not in self.attributes )
		assert( count in ELEMENT_COUNTS )
		assert( not identifier.startswith( '_' ) )
		
		if ':' in identifier:
			print self.identifier
		
		self.attributes[identifier] = count

class TemplateParser( object, sax.handler.ContentHandler ):
	'''Parses the template and builds the parser objects.'''
	__slots__ = ( 'elements', 'elementStack', 'root', 'locator' )
	
	def __init__( self, templateString ):
		'''Creates a new TemplateParser from the templateString.'''
		
		self.elements = {}
		self.elementStack = []
		self.root = None
		self.locator = None
		
		saxparser = sax.make_parser()
		saxparser.setContentHandler( self )
		saxparser.setFeature(sax.handler.feature_namespaces, 1)
		saxparser.parse( StringIO( templateString ) )
		
		self.locator = None
	
	def setDocumentLocator( self, locator ):
		'''SAX Callback: sets the locator so we can get line/column info
		when throwing exceptions.'''
		
		self.locator = locator
	
	def startElementNS( self, name, qname, attributes ):
		'''SAX Callback.'''
		
		if len( self.elementStack ) > 0 and isinstance( self.elementStack[-1], tuple ):
			element, reference = self.elementStack[-1]
			# The last element was a reference: no content permitted
			raise ParseError( E_REF_NOT_EMPTY, self.locator, element, reference )
			
		# Get a Python identifier from the tag name
		identifier = validateIdentifier( name, self.locator )
		
		# Check if this element is a reference
		if attributes.has_key( SIMPLEPARSE_REFERENCE ):
			reference = None
			
			# If the reference is to something that is either non-ascii,
			# or that we don't permit, it is impossible for use to find it
			# in the defined types, so just ignore this error
			try:
				reference = attributes[SIMPLEPARSE_REFERENCE].encode('ascii')
			except UnicodeDecodeError:
				pass
			
			if reference not in self.elements:
				raise ParseError( E_UNDEFINED_REF, self.locator, name, attributes[SIMPLEPARSE_REFERENCE] )
			
			# Validate that there is no content or extra attributes on this item
			for attrName in attributes.keys():
				if attrName[0] != SIMPLEPARSE_NAMESPACE:
					raise ParseError( E_REF_HAS_ATTRIBUTE, self.locator, name, attrName[1] )
			
			# Associate identifier -> reference
			elementInfo = self.elements[reference]
		
		# Not a reference: create a new entry for it
		else:
			if identifier in self.elements:
				raise ParseError( E_ELEMENT_CONFLICT, self.locator, name, None ) 
			assert( identifier not in self.elements )
			
			elementInfo = ElementInfo( name, identifier )
			self.elements[identifier] = elementInfo
		
		if len( self.elementStack ) > 0:
			# Not the initial document element
			count = COUNT_ZERO_OR_ONE
			if attributes.has_key( SIMPLEPARSE_COUNT ):
				count = attributes[SIMPLEPARSE_COUNT]
			
			if count not in ELEMENT_COUNTS:
				# This is not a valid count
				raise ParseError( E_COUNT, self.locator, name, count )
			
			if identifier in self.elementStack[-1].attributes:
				raise ParseError( E_CONFLICT, self.locator, self.elementStack[-1].tagName, identifier )

			self.elementStack[-1].addChild( identifier, elementInfo, count )
		else:
			# This is the initial document element: cannot specify a count
			if attributes.has_key( SIMPLEPARSE_COUNT ):
				raise ParseError( E_INVALID_ATTRIBUTE, self.locator, name, 'count' )
		
		for attrName, value in attributes.items():
			if attrName[0] == SIMPLEPARSE_NAMESPACE:
				if attrName not in SIMPLEPARSE_ATTRIBUTES:
					# This attribute is not permitted in the simpleparse namespace
					raise ParseError( E_INVALID_ATTRIBUTE, self.locator, name, attrName )
				
				# Ignore all attributes in the simpleparse namespace
				continue
			
			# Get a Python identifier from the attribute name
			identifier = validateIdentifier( attrName, self.locator )
			
			count = COUNT_ZERO_OR_ONE
			if value.startswith( 'required' ):
				count = COUNT_ONE
			
			elementInfo.addAttribute( identifier, count )
		
		if attributes.has_key( SIMPLEPARSE_REFERENCE ):
			# Put a tuple on the elementStack, as adding things (text, children) is NOT permitted
			self.elementStack.append( (name, reference) )
		else:
			# Add this element to the stack, as we are going to add things to it
			self.elementStack.append( elementInfo )

	def characters( self, text ):
		'''SAX Callback.'''
		
		if isinstance( self.elementStack[-1], tuple ):
			element, reference = self.elementStack[-1]
			# The last element was a reference: no content permitted
			raise ParseError( E_REF_NOT_EMPTY, self.locator, element, reference )
		
		if not text.isspace():
			# If the template contains non-whitespace text,
			# Enable text on the element
			
			self.elementStack[-1].containsText = True

	def endElementNS( self, name, qname ):
		'''SAX Callback.'''
		
		element = self.elementStack.pop()
		
		# Ignore references
		if isinstance( element, tuple ):
			assert len( self.elementStack ) > 0
		else:
			# After the final endElement (end of document),
			# this will set the root correctly
			assert( element.tagName == name )
			self.root = element
			self.root.createObjectType()

class SimpleXMLParser( object, sax.handler.ContentHandler ):
	'''Parses XML document into Python objects.'''
	
	__slots__ = ( 'root', 'stack', 'locator', )
	
	def setDocumentLocator( self, locator ):
		'''SAX Callback. Stores the locator so we can get line/column
		info when throwing exceptions.'''
		
		self.locator = locator
	
	def startElementNS( self, name, qname, attributes ):
		'''SAX Callback.'''
		
		if len( self.stack ) == 0:
			# The initial doctag
			if name != self.root.tagName:
				raise ParseError( E_INVALID_DOCTAG, self.locator, name, self.root.tagName )
			
			self.root.create()
			self.stack.append( self.root )
			
		else:
			# This is a child element
			if len( self.stack[-1].children ) == 0:
				# This element does not permit children
				raise ParseError( E_NO_CHILDREN, self.locator, self.stack[-1].tagName, name )
			
			# The python identifier is the tag name
			# If the tag is not ascii it will not match and will raise "not found" exceptions
			identifier = name[1]
			
			if identifier not in self.stack[-1].children:
				# This element is not permitted as a child
				raise ParseError( E_INVALID_CHILD, self.locator, self.stack[-1].tagName, name )
			
			# Create an instance of the child element
			childElement = self.stack[-1].children[identifier]
			childElement.create()
	
			# Switch the context to the child element: This enables
			# exception cleanup to work
			last = self.stack[-1]
			self.stack.append( childElement )
			
			# Check if we can have at most one element
			count = last.childrenCount[identifier]
			if count == COUNT_ZERO_OR_ONE or count == COUNT_ONE:
				# Is there already an instance for this child?
				if getattr( last.instance, identifier ) != None:
					line = self.locator.getLineNumber()
					column = self.locator.getColumnNumber()
					raise ParseError( E_MULTIPLE, self.locator, last.tagName, name )
				# Store a reference to the child instance
				setattr( last.instance, identifier, childElement.instance )	
			else:
				assert( count == COUNT_ZERO_OR_MORE or count == COUNT_ONE_OR_MORE )
				# Append a reference to the child instance
				getattr( last.instance, identifier).append( childElement.instance )
		
		# Now set any attributes that are specified
		for attrName, value in attributes.items():
			# The python identifier is the attribute name
			# If it is not ascii it will not match and will raise "not found" exceptions
			identifier = attrName[1]
			
			if identifier not in self.stack[-1].attributes:
				raise ParseError( E_INVALID_ATTRIBUTE, self.locator, name, attrName )
			else:
				assert( getattr( self.stack[-1].instance, identifier ) == None )
				setattr( self.stack[-1].instance, identifier, value )
		
		# Check for any required attributes that are not specified
		for identifier, count in self.stack[-1].attributes.iteritems():
			if count == COUNT_ONE and getattr( self.stack[-1].instance, identifier ) == None:
				raise ParseError( E_MISSING_ATTRIBUTE, self.locator, name, identifier )
	
	def endElementNS( self, name, qname ):
		'''SAX Callback.'''
		
		assert( self.stack[-1].instance != None )
		
		# Verify that we have any required children
		for childName, info in self.stack[-1].children.iteritems():
			count = self.stack[-1].childrenCount[childName]
			
			if (count == COUNT_ONE and getattr( self.stack[-1].instance, childName ) == None) or (count == COUNT_ONE_OR_MORE and len( getattr( self.stack[-1].instance, childName ) ) == 0):
				raise ParseError( E_MISSING, self.locator, self.stack[-1].tagName, info.tagName )
				
		
		# Forget the instance and clear the stack
		if len( self.stack ) > 1:
			# Don't forget the root instance
			self.stack[-1].instance = None
		self.stack.pop()
	
	def characters( self, content ):
		'''SAX Callback.'''
		
		if self.stack[-1].containsText:
			self.stack[-1].instance._text += content
		elif not content.isspace():
			# Non-whitespace causes an error if the document does not contain text
			raise ParseError( E_TEXT, self.locator, self.stack[-1].tagName, None )
	
	def __init__( self, templateString ):
		'''Creates a new parser from the template in templateString. If there
		is any error in the template, this will raise a ParseError exception.'''
		
		# Parse the template
		template = TemplateParser( templateString )
		self.root = template.root
		
		self.stack = []
		self.locator = None
	
	def parse( self, documentString ):
		'''Returns the Python object tree representing tho document in
		documentString. If there is any problem, a ParseError exception
		will be raised.'''
		
		assert( len( self.stack ) == 0 )
		
		saxparser = sax.make_parser()
		saxparser.setContentHandler( self )
		saxparser.setFeature(sax.handler.feature_namespaces, 1)
		
		try:
			saxparser.parse( StringIO( documentString ) )
		except:
			# On exception, drop the stack so we can reuse this parser
			for s in self.stack:
				s.instance = None
			self.stack = []
			
			raise
		
		# The last endElement should empty the stack
		assert( len( self.stack ) == 0 )
		
		# Forget the locator
		self.locator = None
		
		instance = self.root.instance
		self.root.instance = None
		return instance

def printObjectTree( obj, levelString='\n .' ):
	'''Prints a tree representation of obj. Useful for debugging.'''
	
	print 'element', obj.__class__.__name__ ,

	for identifier in obj.__slots__:
		child = getattr( obj, identifier )
		print '%s%s =' %( levelString, identifier ),
		if isinstance( child, ElementBase ):
			printObjectTree( child, levelString + ' .' )
		elif isinstance( child, list ):
			if len(child) == 0:
				print '[]',
			else:
				print '[...]:',
				for i, value in enumerate( child ):
				    #~ print '%s  [%d] of' % ( levelString, i ),
				    printObjectTree( value, levelString + ' [%d].' % i )
		else:
			if child is None:
				print 'None',
			elif identifier == '_text':
			    print '"%s"' % child.encode('unicode_escape') ,
			else:
				s = u'attribute %s "%s" ' % ( identifier, child )
				print 'attribute %s "%s" ' % ( identifier, child ),

# The default main routine that can be useful for debugging
if __name__ == '__main__':
	import sys
	if len( sys.argv ) != 3:
		print 'Usage: simplexmlparse.py <template file> <document file>'
		sys.exit( 1 )
	
	
	print 'Parsing template...'
	f = file(sys.argv[1])
	templateString = f.read()
	f.close()
	
	parser = SimpleXMLParser( templateString )
	
	print 'Parsing XML document...'
	f = file(sys.argv[2])
	docString = f.read()
	f.close()
	
	docObj = parser.parse( docString )
	print 'Python data structure:'
	printObjectTree( docObj )
