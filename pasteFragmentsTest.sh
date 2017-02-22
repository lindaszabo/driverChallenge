#!/bin/bash

# $1 = directory containing subdirectories for each test
# Each subdirectory contains test.fa with the input fasta file and result.txt with the expected output on the first line

# For each subdirectory, calls driver_challenge.py with test.fa,
# compares script output to expected output and prints PASS / FAIL message. 

# usage: sh pasteFragmentsTest.sh tests

PARENT_DIR=$1

for dir in ${PARENT_DIR}/*/
do
  if [[ -d $dir ]]
  then
    echo "python pasteFragments.py -f ${dir}test.fa"
    MY_OUTPUT=`python pasteFragments.py -f ${dir}test.fa`
    EXPECTED_OUTPUT=`head -n 1 ${dir}result.txt`

    if [ "${MY_OUTPUT}" == "${EXPECTED_OUTPUT}" ]
    then
      echo "PASS\n"
    else
      echo "FAIL\n"
    fi
  fi
done
  
