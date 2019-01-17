#!/bin/bash

fname=$1
nmems=$2
members=$3
rtime=$4
nbrhd=$5
proboutfile=$6

echo $nmems > $fname
echo $members >> $fname
echo $rtime >> $fname
echo $nbrhd >> $fname
echo $proboutfile >> $fname

exit 0 
