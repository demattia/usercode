#!/bin/bash

export dirlist="`ls -1d workdirs/*analysis_* `"

# check for missing histogram files
for i in $dirlist ; do
 if [ ! -e $i/histograms.root ] ; then echo missing histograms.root: $i ; fi
done
echo

# check for missing prefilter files
for i in $dirlist ; do 
 if [ ! -e $i/prefilter.root ] ; then echo missing prefilter.root: $i ; fi
done
echo

# check for logfiles that were not compressed
echo stdout files not zipped: `find $dirlist -name "*.o*"|grep -v .gz|wc`
echo stderr files not zipped: `find $dirlist -name "*.e*"|grep -v .gz|wc`
echo

# check for exceptions in general
(for i in `find $dirlist -name "*.e*.gz"` ; do echo $i: `gunzip -c $i|grep -ci exception`|grep -v ": 0" ; done)>EXCEPTIONS
cat EXCEPTIONS|cut -d ":" -f 1 >lgtmp
mv lgtmp EXCEPTIONS

# most frequent individual reason: file open failures on dcache
(for i in `find $dirlist -name "*.e*.gz"` ; do echo $i: `gunzip -c $i|grep -ci "User timeout canceled session"`|grep -v ": 0" ; done)>EXCEPTIONS.DCACHE
cat EXCEPTIONS.DCACHE|cut -d ":" -f 1 >lgtmp
mv lgtmp EXCEPTIONS.DCACHE
echo number of jobs with dcache file open timeouts: `wc EXCEPTIONS.DCACHE`

# other exceptions:
diff EXCEPTIONS EXCEPTIONS.DCACHE|grep "<"|cut -c 3- > EXCEPTIONS.OTHER
echo number of jobs with other exceptions: `wc EXCEPTIONS.OTHER`

# jobs running into timeout, usually because of slow file open attempts
(for i in `find $dirlist -name "*.e*.gz"` ; do echo $i: `gunzip -c $i|grep -ci "job killed"`|grep -v ": 0" ; done)>EXCEPTIONS.JOBKILLED
cat EXCEPTIONS.JOBKILLED|cut -d ":" -f 1 >lgtmp
mv lgtmp EXCEPTIONS.JOBKILLED
echo "number of jobs that were killed (not included in exception count!):" `wc EXCEPTIONS.JOBKILLED`

# jobs running into timeout, usually because of slow file open attempts - but without error message
(for i in `find $dirlist -name "*.e*.gz"` ; do echo $i: `gunzip -c $i|grep -c "TrigReport Events total"`|grep -v ": 1" ; done)>EXCEPTIONS.JOBKILLED2
cat EXCEPTIONS.JOBKILLED2|cut -d ":" -f 1 >lgtmp
mv lgtmp EXCEPTIONS.JOBKILLED2
echo "number of jobs that were killed (not included in exception count!):" `wc EXCEPTIONS.JOBKILLED2`

# jobs that simply crashed
(for i in `find $dirlist -name "*.e*.gz"` ; do echo $i: `gunzip -c $i|grep -c "Segmentation fault"`|grep -v ": 0" ; done)>EXCEPTIONS.SEGFAULT
cat EXCEPTIONS.SEGFAULT|cut -d ":" -f 1 >lgtmp
mv lgtmp EXCEPTIONS.SEGFAULT
echo "number of jobs that had segfaults (not included in exception count!):" `wc EXCEPTIONS.SEGFAULT`

# write script to resubmit jobs with file open errors
echo "#!/bin/bash" > RESUBMIT
echo "export curdir=${PWD}" >> RESUBMIT
for i in `cat EXCEPTIONS`; do
    echo $i
    export jobname=`echo $i|cut -d "/" -f 3|cut -d "e" -f 1|sed -e "s/sh./sh/g"`
    export dirname=workdirs/`echo $i|cut -d "/" -f 2`
    echo rm $i >> RESUBMIT
    echo rm `echo $i|sed -e "s/.sh.e/.sh.o/g"` >> RESUBMIT
    echo cd $dirname >> RESUBMIT
    echo qsub -q prod $jobname >> RESUBMIT
    echo 'cd $curdir' >> RESUBMIT
    echo sleep 0.2 >> RESUBMIT
done

# add jobs to RESUBMIT that for whatever reasons were not even submitted
echo "####### missing jobs follow #######" >> RESUBMIT
for i in $dirlist ; do
    export dirname=workdirs/`echo $i|cut -d "/" -f 2`
    ./HarderAnalysis/DisplacedDileptons/test/find_missing_jobs.py $dirname >> RESUBMIT
done
