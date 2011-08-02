# Split the job and run

import os

eventsPerJob = 10000
totalEvents = 40000000
# eventsPerJob = 1000
# totalEvents = 1000
cmsswDir = '/uscmst1/prod/sw/cms/slc5_amd64_gcc434/cms/cmssw/CMSSW_4_2_3/'

i = 0
while i*eventsPerJob < totalEvents:
    fileNameHLT = 'newFile_HLT_'+str(i)+'.py'
    fileNameReco = 'newFile_Reco_'+str(i)+'.py'
    scriptName = 'Run_'+str(i)+'.csh'
    condorScriptName = 'runOnCondor_'+str(i)
    os.system('cat hltSetupData.py | sed \'s/TOTEVENTS/'+str(eventsPerJob)+'/\' | sed \'s/SKIPEVENTS/'+str(i*eventsPerJob)+'/\' | sed \'s/NUMBER/'+str(i)+'/\' > '+fileNameHLT)
    os.system('cat rereco_pp.py | sed \'s/NUMBER/'+str(i)+'/\' > '+fileNameReco)
    os.system('cat Run.csh | sed \'s-CMSSWDIR-'+cmsswDir+'-\' | sed \'s/HLTFILENAME/'+fileNameHLT+'/\' | sed \'s/RECOFILENAME/'+fileNameReco+'/\' > '+scriptName)
    os.system('chmod +x '+scriptName)
    os.system('cat runOnCondor | sed \'s/SCRIPT/'+scriptName+'/\' | sed \'s/HLTCFG/'+fileNameHLT+'/\' | sed \'s/RECOCFG/'+fileNameReco+'/\' > '+condorScriptName)
    os.system('condor_submit '+condorScriptName)
    # os.system('cd '+cmsswDir+'; cmsenv; cd -; cmsRun '+fileName)
    i += 1

