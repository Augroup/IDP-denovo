#!/usr/bin/env bash
SDIR="./"
MYTDIR="$SDIR/test_data"
MYODIR="$SDIR/test_output"
MYTEMPDIR="$SDIR/test_temp"
echo $MYTDIR
MYCMD="$SDIR/bin/idpdenovo.py $MYTDIR/scaffold_for_test $MYTDIR/lr_for_test $MYTDIR/sr_1.fa $MYTDIR/sr_2.fa -o $MYODIR --specific_tempdir $MYTEMPDIR --threads 24"
echo $MYCMD
$MYCMD
