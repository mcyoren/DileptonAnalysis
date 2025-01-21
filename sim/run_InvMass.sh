#!/bin/bash

Green='\033[0;32m'
NC='\033[0m' # No Color

echo -e "${Green}$1${NC}"

start=`date +%s`

root -l -b << EOF
    gSystem->Load("MyEvent_OLD_C.so")
    gSystem->Load("../AnaTrain/Run14AuAuLeptonComby/MyML_C.so")
    gSystem->Load("../AnaTrain/Run14AuAuLeptonComby/MyEvent_C.so")
    gSystem->Load("InvMass_C.so")
    InvMass("$1",$2,$3,$4)
EOF

end=`date +%s`
duration=`expr $end - $start`
sec=60
mysec=$(($duration % $sec))
mymin=$(($duration / $sec))
myhours=$(($mymin / $sec))

echo -e "${Green}Accumulated time: $myhours hours $mymin minutes $mysec seconds${NC}"