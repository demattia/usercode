# Split the job and run

import os

eventsPerJob = 100000
totalEvents = THEEVENTS
# cmsswDir = '/uscms/home/demattia/TrackingEfficiencyFromCosmics/CMSSW_4_2_2/src/'
cmsswDir = '/uscmst1/prod/sw/cms/slc5_amd64_gcc434/cms/cmssw/CMSSW_4_2_2/'

i = 0
while i*eventsPerJob < totalEvents:
    fileName = 'newFile_'+str(i)+'.py'
    scriptName = 'Run_'+str(i)+'.csh'
    condorScriptName = 'runOnCondor_'+str(i)
    os.system('cat rereco_pp.py | sed \'s/TOTEVENTS/'+str(eventsPerJob)+'/\' | sed \'s/SKIPEVENTS/'+str(i*eventsPerJob)+'/\' | sed \'s/NUMBER/'+str(i)+'/\' > '+fileName)
    os.system('cat Run.csh | sed \'s-CMSSWDIR-'+cmsswDir+'-\' | sed \'s/FILENAME/'+fileName+'/\' > '+scriptName)
    os.system('chmod +x '+scriptName)
    os.system('cat runOnCondor | sed \'s/SCRIPT/'+scriptName+'/\' | sed \'s/CFG/'+fileName+'/\' > '+condorScriptName)
    os.system('condor_submit '+condorScriptName)
    # os.system('cd '+cmsswDir+'; cmsenv; cd -; cmsRun '+fileName)
    i += 1

