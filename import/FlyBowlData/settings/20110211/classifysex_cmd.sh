#!/bin/bash

EXPDIR=$1
SCRIPTDIR=`dirname "$0"`
SCRIPT=$SCRIPTDIR/run_FlyBowlClassifySex.sh 
MCR=/groups/branson/bransonlab/projects/olympiad/MCR/v714
PROTOCOL=20110211

$SCRIPT $MCR $EXPDIR analysis_protocol $PROTOCOL