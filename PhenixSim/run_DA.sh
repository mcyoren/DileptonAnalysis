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
DIR=`printf "%07d" $INPUT`


mkdir $DIR
pushd $DIR
cp /phenix/plhf/mitran/Simul/ETA/pi0_sim/run_DA_sim.C .

export oscar_path=/phenix/plhf/mitran/Simul/ETA/output_pi0/single/$NAME.oscar.particles.dat

if [[ $2 -gt 0 ]]
    then
        echo "embedding"
        for icent in {0..3}
        do  
            echo $icent
            export dst_path=/phenix/plhf/mitran/Simul/ETA/output_pi0/embed/embed_pi0-0-20GeV-$icent'_'$NAME.root
            export out_path=/phenix/plhf/mitran/Simul/ETA/output_pi0/Trees/Tree_$icent'_'$NAME.root
            root -b -q 'run_DA_sim.C('$INPUT',"'$dst_path'","'$out_path'","'$oscar_path'","AnalysisTree")'
        done 
    else 
        echo "no embedding"
        export dst_path=/phenix/plhf/mitran/Simul/ETA/output_pi0/dst/dst_out_pi0_$NAME.root
        export out_path=/phenix/plhf/mitran/Simul/ETA/output_pi0/Trees/Tree$NAME.root
        root -b -q 'run_DA_sim.C('$INPUT',"'$dst_path'","'$out_path'","'$oscar_path'","AnalysisTree0")'
    fi

popd
rm -r $DIR
