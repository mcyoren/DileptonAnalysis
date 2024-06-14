cd pi0_sim

if test -d /phenix/plhf/mitran/Simul/ETA/output_pi0/Trees/; then
echo "dst dir exist"
else 
echo "creating dst dir"
mkdir /phenix/plhf/mitran/Simul/ETA/output_pi0/Trees/
fi 

#cd build/
#make 
#make install
#cd ../
cd /phenix/plhf/mitran/Simul/ETA/output/

shift=400
INPUT=$(( $1 + $shift ))
echo $INPUT
NAME=`printf "%05d" $INPUT`
DIR=`printf "%08d" $INPUT`


mkdir $DIR
pushd $DIR
cp /phenix/plhf/mitran/Simul/ETA/pi0_sim/SimFrac.C .
cp /phenix/plhf/mitran/Simul/ETA/pi0_sim/SimFrac.h .

export tree0_path=/phenix/plhf/mitran/Simul/ETA/output_pi0/Trees/Tree$NAME.root


echo "sim frac"
for icent in {0..3}
do  
    echo $icent
    export tree1_path=/phenix/plhf/mitran/Simul/ETA/output_pi0/Trees/Tree_$icent'_'$NAME.root
    export out_path=/phenix/plhf/mitran/Simul/ETA/output_pi0/Trees/Tree_simfrac_$icent'_'$NAME.root
    root -b -q 'SimFrac.C("'$tree0_path'","'$tree1_path'","'$out_path'")'
done 

popd
rm -r $DIR
