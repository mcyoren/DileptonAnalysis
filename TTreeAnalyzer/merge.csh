#!/bin/bash

start=$((1+${1}*25))
finish=$(($start+24))

  N=`sed -n "$start,$finish p" < output/runs.txt`

echo $N
hadd output/temp/merge_${1}.root $N
