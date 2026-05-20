#!/bin/csh

source /opt/phenix/core/bin/phenix_setup.csh -n

g++ -O2 -o ConvertPythia6PairsToTree \
  ConvertPythia6PairsToTree.cc \
  `root-config --cflags --libs`

./ConvertPythia6PairsToTree $1
