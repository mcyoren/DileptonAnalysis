#!/bin/sh
if [ "$#" -ne 4 ]; then
    echo "You need to run the script with the following 4 parameters: input file, number of threads, number of events, output file name"
    exit 1
fi

echo "input file: $1, number of threads: $2, number of events: $3 output file name: $4"

root -l -b << EOF
    .L ../AnaTrain/Run14AuAuLeptonComby/MyML.C+
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