#!/usr/bin/python

## use os module for running shell commands
## use sys module for getting command line args
import os, sys

## default argument values
jobnumber = 0   ## for condor. corresponds to $(Process).
nevents = 0     ## for Fun4All. 0 means run over all events.
pdgcode = "-321"

## get from command line or condor to override defaults
if len(sys.argv) > 1: jobnumber = int(sys.argv[1])
if len(sys.argv) > 2: pdgcode = sys.argv[2]

hlprefix = "/phenix/hl/tujuba/embedding/prod/dsts/DST_MinBias_run7AuAu_Central_200GeV_pro79-"
data82   = "/phenix/hp/data82/phnxhp01/andrew/emb/"
work     = "/phenix/u/workarea/adare/offline/analysis/runembed/"

## The RD DST list is used here to keep track of run/segment numbers
rdlist = open(work + "jobs/realDST.list", "r").readlines()
runseg = rdlist[jobnumber].split("-")[1] + "-" + rdlist[jobnumber].split("-")[2][0:4]
runnumber = rdlist[jobnumber].split("-")[1][4:]

## Define input and output file names.
RDin        = hlprefix + runseg    + ".root"
MCin        = data82   + "DST/dst" + runseg + "_" + pdgcode + ".root"
ntuple_out  = data82   + "eval/nt" + runseg + "_" + pdgcode + ".root"
DST_out     = "" #data82   + "eval/embdst" + runseg + "_" + pdgcode + ".root"
macroname   = work     + "runembed.C"

command = "root -b -q \'%s(%i, \"%s\", \"%s\", \"%s\", \"%s\")\'"\
                % (macroname, nevents, RDin, MCin, DST_out, ntuple_out)

#print(command)

os.system(command)
