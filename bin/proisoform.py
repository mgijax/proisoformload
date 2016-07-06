#!/usr/local/bin/python

'''
#
# proisoform.py
#
#       See www.informatics.jax.org/wiki/index.php/sw:Proisoformload
#
# Inputs:
#
#       ${OB  the monthly file
#       ${OB  the increment file
#
# Outputs/Re
#
# 	The loader format has the following columns:
#	${INFILE_NAME_VOC}
#
#	1. Term (name)
#	2. Accession ID : PR:xxxx
#	3. current
#	4. Abbreviation (blank)
#	5. Definition
#	6. Comment (blank)
#	7. Synonyms
#	8. Secondary Accession IDs (blank)
#
# 	The annotation loader format has the following columns:
#	${ANNOTINPUTFILE}
#
#       1.  Term ID
#       2.  Marker ID
#       3.  J:
#       4.  Evidence Code
#       5.  Inferred From (blank)
#       6.  Qualifier (blank)
#       7.  Editor : proisoformload
#       8.  Date
#       9.  Notes (blank)
#       10. Logical DB Name of Object (blank)
#       11. Properties (external ref&=&)
#
# History:
#
# 06/29/2016	lec
#       TR12349/TR12345/GPAD/GPI/Noctua
#	converted from Mary's perl script (see TR12345 wts directory)
#
'''

import sys 
import os
import loadlib

obo1FileName = ''
obo1File = ''
obo2FileName = ''
obo2File = ''

# vocabulary formatted file
vocFileName = ''
# vocabulary file pointer
vocFile = ''

# annotation formatted file
annotFileName = ''
# annotation file pointer
annotFile = ''

# load J:
loadjnumber = ''

# load provider
loadprovider = ''

# date for annotation file
loaddate = ''

# xref -> "external ref" for annotation
externalRef = 'external ref&=&'

# symbol, prId, name, synonyms
vocLine = '%s\t%s\tcurrent\t\t%s\t\t%s\t\n'

# prId, mgiId, date
annotLine = '%s\t%s\t%s\tIEA\t\t\t%s\t%s\t\t\t%s\n'

#
# class to store the node information
#
class Node:
    def __init__ (self, prId):
        self.prId = prId
	self.symbol = ''
	self.name = ''
	self.isMouse = 0
	self.synonym = []
	self.parentId = []
	self.mgiId = []
	self.xref = []
		
    def toString (self):
        return '%s|%s|%s|%s|%s|%s\n' % (self.prId, self.symbol, self.parentId, self.mgiId, self.ancestorId, str(self.isMouse))

#
# the nodes that we will actually use
# to determine the prId->mgiId association
#
# if isMouse = 1, then this node will be included in the vocabulary/annotation files.
#
# if isMouse = 0, then the node is needed to navigate thru the hierarchy
# but will not be included in the vocabulary/annotation files.
#
nodeLookup = {}

#
# due to the excluded nodes (see skipPR), some nodes will not exist in nodeLookup.
# however, we still need the prId->mgiId associations
#
mgiLookup = {}

#
# Purpose: Initialization
#
def initialize():

    global obo1FileName, obo1File
    global obo2FileName, obo2File
    global vocFileName, vocFile
    global annotFileName, annotFile
    global loadjnumber, loadprovider, loaddate

    #
    # open files
    #

    obo1FileName = os.environ['OBO1FILE']
    obo2FileName = os.environ['OBO2FILE']
    vocFileName = os.environ['INFILE_NAME_VOC']
    annotFileName = os.environ['ANNOTINPUTFILE']

    loadjnumber = os.environ['JNUMBER']
    loadprovider = os.path.basename(os.environ['PROISOFORMLOAD'])
    loaddate = loadlib.loaddate

    obo1File = open(obo1FileName, 'r')
    obo2File = open(obo2FileName, 'r')
    vocFile = open(vocFileName, 'w')
    annotFile = open(annotFileName, 'w')

    return

#
# Purpose: Read/Process obo files and generate Vocabulary/Annotation files
#
def processOBO(oboFile):

    startTag = '[Term]'
    endTag = ''
    prTag = 'id: PR:'
    idTag = 'id: '
    skipPR = ['PR:000029032','PR:000000001','PR:000018263']

    startTerm = 0
    foundTerm = 0
    foundRelationship = 0
    addToLookup = 1

    #
    # load terms into nodeLookup
    #

    for line in oboFile.readlines():

	line = line.strip()

	# 
	# start of term
	#
	if line.find(startTag) == 0:
	    startTerm = 1
	    foundTerm = 0
	    foundRelationship = 0
            addToLookup = 1
	    continue

	#
	# end of term
	#
	if startTerm and len(line) == 0:

	    #
	    # if nodeLookup already exist for this term, then use it
	    #
	    try:
		if foundTerm and foundRelationship and addToLookup:
	            nodeLookup[n.prId] = n
	    except:
	    	pass

	    startTerm = 0
	    foundTerm = 0
	    foundRelationship = 0
            addToLookup = 1
	    continue

	if not startTerm:
		continue

 	if line.find(prTag) == 0:

	    ignoreit, prId = line.split(idTag)

	    #
	    # prId may be in both obo files
	    #
	    if prId in nodeLookup:
	        if prId in ('PR:A2CG49-2'):
		    p = nodeLookup[prId]
	        n = nodeLookup[prId]
	    else:
	        n = Node(prId)
                n.prId = prId

	    # if prId like 'xxxx-1', then attach 'xxxx' as parent
	    try:
	    	nodeParent, ignoreit = prId.split('-')
		if nodeParent not in n.parentId:
    	            n.parentId.append(nodeParent)
	    except:
	        pass

	    foundTerm = 1
	    continue

	#
	# could not find Term
	#
	if not foundTerm:
	    continue

	#
	# found proper term, continue
	#

	if line.find('name:') == 0:
	    n.name = line[6:]

	#
	# use this as the symbol 
        # synonym: "xxxx" EXACT PRO-short-label [PRO:DNx]
	#
	elif line.find('synonym:') == 0 and line.find('EXACT PRO-short-label') >= 0:
            tokens = line.split(' ')
            n.symbol = tokens[1].replace('"','')
	    n.synonym.append(tokens[1].replace('"',''))

	#
        # synonym: "xxxxx" RELATED []
        # synonym: "xxxx" EXACT []
	#
	elif line.find('synonym:') == 0:
	    if line.find('EXACT [') >= 0 or line.find('RELATED [') >= 0:
                tokens = line.split(' ')
	        n.synonym.append(tokens[1].replace('"',''))

        #
        # foundRelationship
	#
	# list of typs of "tags" that need to be included in nodeLookup
	#
	# is_a
	# intersection_of: PR:
	# intersection_of: derives_from
	# relationship: derives_from
	# comment: dervies_from
	# intersection_of: has_gene_template MGI:
        # relationship: has_gene_template MGI:
        #

	# is_a
	elif line.find('is_a: PR:') == 0:
	    tokens = line.split(' ')
	    parentId = tokens[1]

	    if parentId not in skipPR:
		if parentId not in n.parentId:
    	            n.parentId.append(parentId)
		foundRelationship = 1

	# intersection_of: PR:
	elif line.find('intersection_of: PR:') == 0:
	    tokens = line.split(' ')
	    parentId = tokens[1]

	    if parentId not in skipPR:
		if parentId not in n.parentId:
    	            n.parentId.append(parentId)
		foundRelationship = 1

	# intersection_of: derives_from
	# relationship: derives_from
	elif line.find('intersection_of: derives_from PR:') == 0 or line.find('relationship: derives_from PR:') == 0:

	    tokens = line.split(' ')
	    parentId = tokens[2]

	    if parentId not in skipPR:
		if parentId not in n.parentId:
    	            n.parentId.append(parentId)
		foundRelationship = 1

	# comment: dervies_from
	elif line.find('derives_from PR:') >= 0:
	    tokens1 = line.split(' ')
	    for t in tokens1:
	    	if t.find('PR:') == 0:
			tokens2 = t.split('.')

	    parentId = tokens2[0]
	    if parentId not in skipPR:
		if parentId not in n.parentId:
    	            n.parentId.append(parentId)
		foundRelationship = 1

	#
	# to find MGI:xxxx
	# intersection_of: has_gene_template MGI:
        # relationship: has_gene_template MGI:
	#
	elif line.find('intersection_of: has_gene_template MGI:') == 0 \
	     or line.find('relationship: has_gene_template MGI:') == 0:

	    tokens = line.split(' ')
	    mgiId = tokens[2]
	    if mgiId not in n.mgiId:
    	        n.mgiId.append(mgiId)

	    #
	    # for all parents of this node
	    #    create node of parent (if it does not already exist) 
	    #	 attach mgiId
	    #
	    for p in n.parentId:
	        if p not in nodeLookup:
	            newp = Node(p)
                    newp.prId = p
	        else:
	            newp = nodeLookup[p]
		if mgiId not in newp.mgiId:
	            newp.mgiId.append(mgiId)

	    foundRelationship = 1

        #
        # end foundRelationship
        #

	#
	# xref: UniProtKB:xxx
	#
	elif line.find('xref: UniProtKB:') == 0:
	    tokens = line.split(' ')
	    uniprotId = tokens[1]
	    n.xref.append(uniprotId)

	#
	# isMouse?
	#
	elif line.find('only_in_taxon NCBITaxon:10090') >- 0:
	    n.isMouse = 1

	#
	# if 'only_in_taxon' exists and a non-mouse organism is defined,
	#	then addToLookup = 0
	# else
	#	then addToLookup = 1
	#
	# so this will pick up terms that do not have 'only_in_taxon'
	# and are non-mouse as well as "real" mouse terms
	# this is needed for the finding the correct
	# parent for the MGI id
	#
	elif line.find('only_in_taxon NCBITaxon:') >= 0 and line.find('only_in_taxon NCBITaxon:10090') < 0:
	    addToLookup = 0

	#
	# if obsolete, do not add to lookup
	#
	elif line.find('is_obsolete: true') >= 0:
	    addToLookup = 0

    # last record
    try:
	if foundRelationship and addToLookup:
	    nodeLookup[n.prId] = n
    except:
        pass

    return

#
# Purpose: Recursive function that iterates thru parentId
#          until it finds a parent that contains an mgiId
#
def findMgiIdByParent(n):

    #if n.prId in ('PR:000025728'):
      #print n.toString()

    mgiId = ''

    # for each parentId in the node

    for p1 in n.parentId:

	# get the node of the parentId
	# the parent may not be in nodeLookup
	# so check mgiLookup...
	if p1 not in nodeLookup:
	    if p1 in mgiLookup:
	        mgiId = mgiLookup[p1]
		return mgiId
	    else:
	        continue

	p2 = nodeLookup[p1]

	# does the node contain an mgiId?
	# if not, keep looking

	if len(p2.mgiId) == 0:
	    mgiId = findMgiIdByParent(p2)

	# else, we are done
	else:
	    mgiId = p2.mgiId
            return mgiId

    return mgiId

#
# Purpose: To print the nodeLookup information into the correct output files
#
def printFiles():

    #
    # build mgiLookup (see comments above)
    # 
    for r in nodeLookup:

	n = nodeLookup[r]

	for p in n.parentId:
	    if len(n.mgiId) > 0:
	        if p not in mgiLookup:
	            mgiLookup[p] = []
		    mgiLookup[p].append(n.mgiId[0])

    for r in nodeLookup:

	n = nodeLookup[r]

	#if n.prId in ('PR:000037266', 'PR:Q62507'):
	    #print genericLookup['PR:000002187']
    	    #print n.toString()

	# find mgiId of term
	# uses 'findMgiIdByParent() to iterate thru each parentId

	if len(n.mgiId) == 0:
	    mgiId = findMgiIdByParent(n)
        else:
	    mgiId = n.mgiId

	#
	# if no mgiId was found, skip it
	#
	if len(mgiId) == 0 or not n.isMouse:
	    continue

	#
	# if no "symbol" given, then use "name"
	#
	if len(n.symbol) == 0:
	    symbol = n.name
        else:
	    symbol = n.symbol

	#
	# synonyms
	#
	synonyms = '|'.join(n.synonym)

	vocFile.write(vocLine % (symbol, n.prId, n.name, synonyms))

	#
	# property
	#
	property = ''
	if len(n.xref) > 0:
	    property = externalRef + n.xref[0]

        annotFile.write(annotLine % (n.prId, mgiId[0], loadjnumber, loadprovider, loaddate, property))

    return

#
# Purpose: Close files
#
def closeFiles():

    obo1File.close()
    obo2File.close()
    vocFile.close()
    annotFile.close()

    return

#
# main
#

initialize()
processOBO(obo1File)
processOBO(obo2File)
printFiles()
closeFiles()

