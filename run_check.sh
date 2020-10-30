#!/usr/bin/env bash

pyver=$(python -c "import sys; print('{0[0]}.{0[1]}'.format(sys.version_info))")
echo "python version is "$pyver
if [ "$pyver" != "2.7" ]; then
  echo "ERROR $pyver is not python 2.7"
  exit 64
fi


SDIR="./"
MYTDIR="$SDIR/test_data"
echo $MYTDIR
MYCMD="$SDIR/bin/idpdenovo.py $MYTDIR/scaffold_for_test $MYTDIR/lr_for_test $MYTDIR/sr_1.fa $MYTDIR/sr_2.fa --test -o TESTING"
echo $MYCMD
$MYCMD
