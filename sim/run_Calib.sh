#!/bin/sh

echo "$1"

start=`date +%s`

root -l -b << EOF
    gSystem->Load("../AnaTrain/Run14AuAuLeptonComby/MyML_C.so")
    gSystem->Load("../AnaTrain/Run14AuAuLeptonComby/MyEvent_C.so")
    gSystem->Load("Calib_C.so")
    Calib("$1",$2,$3,$4)
EOF

end=`date +%s`
duration=`expr $end - $start`
sec=60
mysec=$(($duration % $sec))
mymin=$(($duration / $sec))
myhours=$(($mymin / $sec))

echo "Accumulated time:" $myhours "hours" $mymin "minutes" $mysec "seconds"