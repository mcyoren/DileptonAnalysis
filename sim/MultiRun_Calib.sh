#!/bin/sh

echo "$1 $2 $3"

root -l -b << EOF
    .L ../AnaTrain/Run14AuAuLeptonComby/MyEvent.C+
    .L Calib.C+
EOF

for i in `seq 1 $2`
do
    ./run_Calib.sh $1 $i $2 &
done
wait
hadd -f -k -O $3 my-kek*
rm my-kek*
cp $3 output/my-kek.root
echo "done!"