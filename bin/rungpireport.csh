#!/bin/csh -f

#
# rungpireport.csh
#
# Script to run GPI report from public report product (reports_db)
#

cd `dirname $0` && source ${PUBRPTS}/Configuration

cd ${PUBRPTS}/daily
./GO_gpi.py
