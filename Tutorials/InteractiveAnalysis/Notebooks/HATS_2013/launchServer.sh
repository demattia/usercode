source /home/ipython1/soft/setup.sh
export LD_LIBRARY_PATH=/afs/cern.ch/cms/slc5_amd64_gcc462/cms/cmssw/CMSSW_5_3_8/lib/slc5_amd64_gcc462:/afs/cern.ch/cms/slc5_amd64_gcc462/cms/cmssw/CMSSW_5_3_8/external/slc5_amd64_gcc462/lib:$LD_LIBRARY_PATH
export CMSSW_BASE=/afs/cern.ch/cms/slc5_amd64_gcc462/cms/cmssw/CMSSW_5_3_8
export CMSSW_RELEASE_BASE=/afs/cern.ch/cms/slc5_amd64_gcc462/cms/cmssw/CMSSW_5_3_8
export PYTHONPATH=/afs/cern.ch/cms/slc5_amd64_gcc462/cms/cmssw/CMSSW_5_3_8/python:/afs/cern.ch/cms/slc5_amd64_gcc462/cms/cmssw/CMSSW_5_3_8/lib/slc5_amd64_gcc462:$PYTHONPATH
nohup ipython notebook --profile=myserver &