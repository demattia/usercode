#!/usr/bin/env python

import os

def submit(leptype):
    counter = 0
    os.system('mkdir -p '+leptype)
    for mass in range(20,505,5):
        # if counter == 2:
        #     break
        counter += 1
        jobFileName = leptype+'/'+leptype+'_job_'+str(mass)+'.sh'
        jobFile = open(jobFileName, 'w')
        jobFile.write('source '+os.getcwd()+'/cmslpc_standalone_setup.sh\n')
        jobFile.write('cd '+os.getcwd()+'\n')
        jobFile.write('python main.py '+str(mass)+' '+leptype+'\n')
        jobFile.close()
        os.system('chmod +x '+jobFileName)

        condorFileName = leptype+'/runOnCondor_'+str(mass)
        condorFile = open(condorFileName, 'w')
        condorFile.write('universe = vanilla\n')
        condorFile.write('Executable = '+jobFileName+'\n')
        condorFile.write('Requirements = Memory >= 199 &&OpSys == "LINUX"&& (Arch != "DUMMY" )&& Disk > 1000000\n')
        condorFile.write('Should_Transfer_Files = YES\n')
        condorFile.write('Transfer_Input_Files = '+jobFileName+'\n')
        condorFile.write('WhenToTransferOutput = ON_EXIT\n')
        # condorFile.write('Log = condor_job1_cfg.py_$(Cluster)_$(Process).log\n')
        # condorFile.write('Output = condor_rereco_pp_$(Cluster)_$(Process).stdout\n')
        # condorFile.write('Error = condor_rereco_pp_$(Cluster)_$(Process).stderr\n')
        condorFile.write('notify_user = '+os.environ["USER"]+'@FNAL.GOV\n')
        condorFile.write('Queue 1\n')
        condorFile.close()

        print 'Submitting '+leptype
        os.system('condor_submit '+condorFileName)

submit("Electrons")

submit("Muons")
