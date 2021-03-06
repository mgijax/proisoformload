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
# Columns:
# 
#    name                   required? cardinality   GAF column #  Note
#    DB                     required  1             1
#    DB_Object_ID           required  1             2/17
#    DB_Object_Symbol       required  1             3             Based on the PRO-short-label
#    DB_Object_Name         optional  0 or greater  10
#    DB_Object_Synonym(s)   optional  0 or greater  11            Only exact synonyms are given
#    DB_Object_Type         required  1             12            Will always be 'protein' or 'protein_complex' for PRO terms
#    Taxon                  required  1             13
#    Parent_Object_ID       optional  0 or 1        -             Blank if DB_Object_ID refers to a canonical entity, otherwise it will be an ancestral identifier (in the ontological sense) that refers to the canonical entity class.
#    DB_Xref(s)             optional  0 or greater  -             The gene(s) encoding the DB_Object_ID; multiple genes are separated by pipes
#    Properties             optional  0 or greater  -   
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
# 08/12/2020    lec
#       TR13272/added OUTPUT_GPI2 file
#
# 05/04/2020	lec
#       TR13299/switch to use the PRO/GPI
#
# 06/29/2016	lec
#       TR12349/TR12345/GPAD/GPI/Noctua
#	converted from Mary's perl script (see TR12345 wts directory)
#
'''

import sys 
import os
import db
import mgi_utils
import loadlib

diagFile = ''
curFile = ''
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

# output files
gpi2FileName = ''
gpi2File = ''

# load J:
loadjnumber = ''

# load provider
loadprovider = ''

# date for annotation file
loaddate = ''

# xref -> "external ref" for annotation
externalRef = 'external ref&=&'

# symbol, prId, synonym
vocLine = '%s\t%s\tcurrent\t\t%s\t\t%s\t\n'

# prId, mgiId, date
annotLine = '%s\t%s\t%s\tNAS\t\t\t%s\t%s\t\t\t%s\n'

#
# Purpose: Initialization
#
def initialize():

    global diagFile, curFile
    global gpiFileName, gpiFile
    global vocFileName, vocFile
    global annotFileName, annotFile
    global loadjnumber, loadprovider, loaddate
    global gpi2FileName, gpi2File

    #
    # open files
    #

    diagFileName = os.environ['LOG_DIAG']
    curFileName = os.environ['LOG_CUR']
    gpiFileName = os.environ['GPIFILE']
    vocFileName = os.environ['INFILE_NAME_VOC']
    annotFileName = os.environ['ANNOTINPUTFILE']
    gpi2FileName = os.environ['OUTPUT_GPI2']

    loadjnumber = os.environ['JNUMBER']
    loadprovider = os.path.basename(os.environ['PROISOFORMLOAD'])
    loaddate = loadlib.loaddate

    try:
        diagFile = open(diagFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % diagFileName)
            
    try:
        curFile = open(curFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % curFileName)
            
    try:
        gpiFile = open(gpiFileName, 'r')
    except:
        exit(1, 'Could not open file %s\n' % gpiFileName)
            
    try:
        vocFile = open(vocFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % vocFileName)
            
    try:
        annotFile = open(annotFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % annotFileName)
            
    try:
        gpi2File = open(gpi2FileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % gpi2FileName)
            
    # Log all SQL 
    db.set_sqlLogFunction(db.sqlLogAll)

    # Set Log File Descriptor
    db.set_commandLogFile(diagFileName)

    # Set Log File Descriptor
    diagFile.write('Start Date/Time: %s\n' % (mgi_utils.date()))
    diagFile.write('Server: %s\n' % (db.get_sqlServer()))
    diagFile.write('Database: %s\n' % (db.get_sqlDatabase()))
    diagFile.write('Input File: %s\n' % (gpiFileName))
    curFile.write('Start Date/Time: %s\n\n' % (mgi_utils.date()))

    return

#
# Purpose: Read/Process files and generate Vocabulary/Annotation files
#
def processGPI():

    for line in gpiFile.readlines():

        if line.find('!') >= 0:
            continue

        line = line.strip()
        tokens = line.split('\t')
        prId = tokens[0] + ':' + tokens[1]
        uniprotId = tokens[1]
        symbol = tokens[2]
        name = tokens[3]
        synonym = tokens[4]
        prtype = tokens[5]
        taxon = tokens[6]
        try:
            mgiId = tokens[8].replace('MGI:MGI:', 'MGI:')
        except:
            mgiId = ''
        #
        # mouse only
        #
        if taxon != 'taxon:10090':
            continue

        #
        # missing MGI:xxxx && protein_complex
        #
        if mgiId.find('MGI:') < 0 and prtype == 'protein_complex':
            taxon = taxon.replace('taxon', 'NCBITaxon')
            gpi2File.write(prId + '\t')
            gpi2File.write(symbol + '\t')
            gpi2File.write(name + '\t')
            gpi2File.write(synonym + '\t')
            gpi2File.write('GO:0032991' + '\t')
            gpi2File.write(taxon + '\t')
            gpi2File.write('\t')
            gpi2File.write('\t')
            gpi2File.write('\t')
            gpi2File.write('\t')
            gpi2File.write('\n')

        #
        # missing MGI:xxxx
        #
        if mgiId.find('MGI:') < 0:
            curFile.write('missing MGI:xxxx\n')
            curFile.write(line + "\n\n")
            continue

        #
        # 1 mouse only
        #
        if mgiId.find('|') > 1:
            curFile.write('> 1 mouse MGI:xxx\n')
            curFile.write(line + "\n\n")
            continue

        #
        # check if UniProtKB sequence exists
        #
        results = '''select * from ACC_Accession where _mgitype_key = 19 and accid = '%s' ''' % (uniprotId)
        if len(results) == 0:
            curFile.write('uniprotKB not found in MGI\n')
            curFile.write(line + "\n\n")
            continue

        property = ''
        if prId[3] != '0':
            property = externalRef + 'UniProtKB:' + uniprotId

        vocFile.write(vocLine % (symbol, prId, name, synonym))
        annotFile.write(annotLine % (prId, mgiId, loadjnumber, loadprovider, loaddate, property))

    return

#
# Purpose: Close files
#
def closeFiles():

    diagFile.flush()
    curFile.flush()
    diagFile.write('\n\nEnd Date/Time: %s\n' % (mgi_utils.date()))
    diagFile.close()
    curFile.close()
    gpiFile.close()
    vocFile.close()
    annotFile.close()
    gpi2File.close()

    return

#
# main
#

initialize()
processGPI()
closeFiles()

