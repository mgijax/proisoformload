#!/usr/local/bin/python

'''
#
# proisoform.py
#
#	build vocabulary & annotation output text files by using the PR GPI input files.
#	this data will be used to generate a GO/GPI file (see reports_db/daily/GO_gpi.py)
#
#	PR:xxxx are stored as a simple vocabulary (no DAG)
#	PR:xxxx/MGI:xxxx associations are stored as Annotations
#
#	PR:xxxx has one MGI:xxxx
#	MGI:xxxx can have >= 1 PR:xxxx
#
#       See www.informatics.jax.org/wiki/index.php/sw:Proisoformload
#
# Inputs:
#
#       ${GPIFILE}
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
# 05/04/2020	lec
#       TR13299/swtich to use the PRO/GPI
#
# 06/29/2016	lec
#       TR12349/TR12345/GPAD/GPI/Noctua
#	converted from Mary's perl script (see TR12345 wts directory)
#
'''

import sys 
import os
import loadlib

gpiFileName = ''
gpiFile = ''

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

# symbol, prId
vocLine = '%s\t%s\tcurrent\t\t\t\t\t\n'

# prId, mgiId, date
annotLine = '%s\t%s\t%s\tIEA\t\t\t%s\t%s\t\t\t%s\n'

#
# Purpose: Initialization
#
def initialize():

    global gpiFileName, gpiFile
    global vocFileName, vocFile
    global annotFileName, annotFile
    global loadjnumber, loadprovider, loaddate

    #
    # open files
    #

    gpiFileName = os.environ['GPIFILE']
    vocFileName = os.environ['INFILE_NAME_VOC']
    annotFileName = os.environ['ANNOTINPUTFILE']

    loadjnumber = os.environ['JNUMBER']
    loadprovider = os.path.basename(os.environ['PROISOFORMLOAD'])
    loaddate = loadlib.loaddate

    gpiFile = open(gpiFileName, 'r')
    vocFile = open(vocFileName, 'w')
    annotFile = open(annotFileName, 'w')

    return

#
# Purpose: Read/Process files and generate Vocabulary/Annotation files
#
def processGPI():

    for line in gpiFile.readlines():

        if line[0] == '!':
            continue

        line = line.strip()
        tokens = line[:-1].split('\t')
        prId = tokens[0] + ':' + tokens[1]
        symbol = tokens[2]
        taxon = tokens[6]
        try:
            mgiId = tokens[8].replace('MGI:MGI:', 'MGI:')
        except:
            mgiId = ''

        #
        # mouse only
        #
        if taxon == 'taxon:10090' and mgiId.find('MGI:') >= 0 :
            mgiId = tokens[8].replace('MGI:MGI:', 'MGI:')
            vocFile.write(vocLine % (symbol, prId))
            annotFile.write(annotLine % (prId, mgiId, loadjnumber, loadprovider, loaddate, externalRef + prId.replace('PR:', 'UniProtKB:')))

    return

#
# Purpose: Close files
#
def closeFiles():

    gpiFile.close()
    vocFile.close()
    annotFile.close()

    return

#
# main
#

initialize()
processGPI()
closeFiles()
