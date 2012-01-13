#!/usr/bin/env python

import os, sys

if len(sys.argv) != 2 or (sys.argv[1] != "Electrons" and sys.argv[1] != "Muons"):
    print "Usage: runEventsLimitsJobs Electrons|Muons"
    sys.exit(1)
leptype = sys.argv[1]

def run(leptype):
    counter = 0
    os.system('mkdir -p '+leptype)
    for mass in range(20,505,5):
        # if counter == 2:
        #     break
        counter += 1
        scriptName = leptype+'_job_'+str(mass)+'.sh'
        jobFileName = leptype+'/'+scriptName
        jobFile = open(jobFileName, 'w')
        jobFile.write('source '+os.getcwd()+'/cmslpc_standalone_setup.sh\n')
        jobFile.write('cd '+os.getcwd()+'\n')
        jobFile.write('python main.py '+str(mass)+' '+leptype+'\n')
        jobFile.close()
        os.system('chmod +x '+jobFileName)

        print 'Running '+leptype
        os.system('cd '+leptype+'; ./'+scriptName)

run(leptype)
