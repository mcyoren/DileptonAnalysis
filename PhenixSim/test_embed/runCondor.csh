#!/bin/tcsh

#-----------------------------------------------------------------------#
# set directory paths
#-----------------------------------------------------------------------#    
set runScript=/gpfs/mnt/gpfs02/phenix/hhj/hhj1/ahodges/run14HadEff/scripts/pisaToPisaDST.csh
#echo "runScript= "$runScript

#set condorDir=/gpfs/mnt/gpfs02/phenix/hhj/hhj1/cwong14/run15hadronEff/condor
set condorDir=/gpfs/mnt/gpfs02/phenix/hhj/hhj1/ahodges/run14HadEff/condor
set logDir=$condorDir/log
set jobDir=$condorDir/job

#-----------------------------------------------------------------------#
# set parameters
#-----------------------------------------------------------------------#
set totalEvent=1000000
set nEvent=5000

#set species=piminus
#set processID=0

foreach species (piminus piplus kminus kplus antiproton proton)

set processID=0
while ( (( $processID * $nEvent )) < $totalEvent )

#-----------------------------------------------------------------------#
# set file names
#-----------------------------------------------------------------------#
set fileTag=$species"_"$processID
set jobFile=$jobDir/$fileTag".job"
set logFile=/tmp/ahodges/logs/$fileTag".log"
set errFile=$logDir/$fileTag".err"
set outFile=$logDir/$fileTag".out"
	
if ( -f $jobFile ) then
    rm -f $jobFile
endif

if ( -f $logFile ) then
    rm -f $logFile
endif

if ( -f $errFile ) then
    rm -f $errFile
endif

if ( -f $outFile) then
    rm -f $outFile
endif
#-----------------------------------------------------------------------#
# write condor job scripts
#-----------------------------------------------------------------------#
echo "Universe = vanilla" >> $jobFile 
echo "Notification = Complete" >> $jobFile
echo "Executable = "$runScript >> $jobFile
echo "Arguments = "$processID" "$nEvent" "$species >>$jobFile
echo "Requirements = CPU_Speed >= 1" >> $jobFile
echo "Requirements = CPU_Experiment == "\""phenix"\" >> $jobFile
echo "Rank = CPU_Speed" >> $jobFile
echo "Image_Size = 500000" >> $jobFile
echo "Priority = +20" >> $jobFile
echo "GetEnv = True" >> $jobFile
#echo "Initialdir = /gpfs/mnt/gpfs02/phenix/hhj/hhj1/cwong14/run15hadronEff/" >> $jobFile
echo "Initialdir = /gpfs/mnt/gpfs02/phenix/hhj/hhj1/ahodges/run14HadEff" >> $jobFile
echo "Input = /dev/null" >> $jobFile
echo "Output = "$outFile >> $jobFile
echo "Log = "$logFile >> $jobFile
echo "Error = "$errFile >> $jobFile
##print condor "+Experiment = \"general\"\n";
echo "+Job_Type = "\""cas"\" >> $jobFile
echo "Queue" >> $jobFile

#$counter = $counter +1;

#-----------------------------------------------------------------------#
# submit job
#-----------------------------------------------------------------------#
condor_submit $jobFile

@ processID++

end
end
