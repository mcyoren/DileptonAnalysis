#!/bin/bash

rm /phenix/plhf/mitran/Analysis/Run14AuAuDiLeptonAnalysis/TTreeAnalyzer/FirstIter/output/temp/*

ls /phenix/plhf/mitran/Analysis/Run14AuAuDiLeptonAnalysis/TTreeAnalyzer/FirstIter/output/All/my-* > output/runs.txt

condor_submit run_merge.job

