#!/bin/sh

echo "$1 $2 $3 $4"

root -l -b << EOF
    .L ../AnaTrain/Run14AuAuLeptonComby/MyEvent.C+
    .L Calib.C+
EOF

for i in `seq 1 $2`
do
    ./run_Calib.sh $1 $i $2 $3 &
done
wait
hadd -f -k -O $4 my-kek*
rm my-kek*
cp $4 output/my-kek.root
echo "done!"