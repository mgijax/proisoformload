#!/bin/sh
#
#  proisoformload.sh
###########################################################################
#
#  Purpose:
# 	This script creates Proteoform vocload and assocload
#       input file and invokes vocload and assocload
#
#  Configuration/Inputs/Outputs:
#
#	proisoform.config
#	vocload.config
#	annotload.config
#
#  Exit Codes:
#
#      0:  Successful completion
#      1:  Fatal error occurred
#      2:  Non-fatal error occurred
#

cd `dirname $0`

COMMON_CONFIG=${PROISOFORMLOAD}/proisoform.config

USAGE="Usage: proisoform.sh"

#
# Make sure the common configuration file exists and source it.
#
if [ -f ${COMMON_CONFIG} ]
then
    . ${COMMON_CONFIG}
else
    echo "Missing configuration file: ${COMMON_CONFIG}"
    exit 1
fi

#
# Initialize the log file.
#
LOG=${LOG_DIAG}
rm -rf ${LOG}
touch ${LOG}

#
# Source the DLA library functions.
#
if [ "${DLAJOBSTREAMFUNC}" != "" ]
then
    if [ -r ${DLAJOBSTREAMFUNC} ]
    then
        . ${DLAJOBSTREAMFUNC}
    else
        echo "Cannot source DLA functions script: ${DLAJOBSTREAMFUNC}" | tee -a ${LOG}
        exit 1
    fi
else
    echo "Environment variable DLAJOBSTREAMFUNC has not been defined." | tee -a ${LOG}
    exit 1
fi

#####################################
#
# Main
#
#####################################

#
# createArchive including OUTPUTDIR, startLog, getConfigEnv
# sets "JOBKEY"
preload ${OUTPUTDIR}

echo "Copy GPI files & running perl script to generate vocload & annotload txt files" >> ${LOG_DIAG}
cp -r ${GPIFILE_DOWNLOAD} ${INPUTDIR}

echo "Generating Proteoform input files" >> ${LOG_DIAG}
${PYTHON} ${PROISOFORMLOAD}/bin/proisoform.py &>> ${LOG}
STAT=$?
checkStatus ${STAT} "${PROISOFORMLOAD}/bin/proisform.py"

#
# run annotation load : delete current data
# else vocabulary cannot be deleted (referencial integrity)
#
cd ${OUTPUTDIR}
echo "Running Proteoform annotation deletion load" >> ${LOG_DIAG}
${ANNOTLOAD}/annotload.csh ${PROISOFORMLOAD}/annotloaddelete.config &>> ${LOG}
STAT=$?
checkStatus ${STAT} "${ANNOTLOAD} ${PROISOFORMLOAD}/annotloaddelete.config"

#
# run the vocabulary load
#
cd ${OUTPUTDIR}
echo "Running vocload to load Proteoform" >> ${LOG_DIAG}
${VOCLOAD}/runSimpleIncLoad.sh ${PROISOFORMLOAD}/vocload.config &>> ${LOG}
STAT=$?
checkStatus ${STAT} "${VOCLOAD}/runSimpleIncLoad.sh ${PROISOFORMLOAD}/vocload.config"

#
# run annotation load
#
cd ${OUTPUTDIR}
echo "Running Proteoform annotation load" >> ${LOG_DIAG}
${ANNOTLOAD}/annotload.csh ${PROISOFORMLOAD}/annotload.config &>> ${LOG}
STAT=$?
checkStatus ${STAT} "${ANNOTLOAD} ${PROISOFORMLOAD}/annotload.config"

#
# run postload cleanup and email logs
#
shutDown
