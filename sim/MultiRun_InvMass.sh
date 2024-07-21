#!/bin/sh

echo "$1 $2 $3"

root -l -b << EOF
    .L MyEvent_OLD.C+
    .L ../AnaTrain/Run14AuAuLeptonComby/MyEvent.C+
    .L InvMass.C+
EOF

for i in `seq 1 $2`
do
    ./run_InvMass.sh $1 $i $2 &
done
wait
hadd -f -k -O output/invmass_embed/$3.root my-kek*
rm my-kek*
echo "done!"