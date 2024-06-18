#!/usr/local/bin/tcsh -f

#set list = `cat run421716.txt`
set list = `cat run409471_20.txt`

@ i = 0

foreach file ( $list )

  set runnum = `basename $file .PRDFF| awk -F'-' '{printf "%s", $2}'`
  set seqnum = `basename $file .PRDFF| awk -F'-' '{printf "%s", $3}'`
  set jobno  = "merge_${runnum}_${seqnum}"
#  echo $runnum $seqnum $jobno

  echo "./jobsubmit.csh reco_rawhit.csh $jobno 0 $file"
        ./jobsubmit.csh reco_rawhit.csh $jobno 0 $file

 @ i ++

#  if($i > 0 ) then
#    exit 0;
#  endif
end
