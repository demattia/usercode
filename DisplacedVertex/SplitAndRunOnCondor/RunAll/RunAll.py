import os
import re
                    
print "Submitting all jobs"

runList = [163852, 163851, 163850, 163847, 162911, 158270, 158197, 158097, 158095, 157989, 157957, 157956, 157953]
eventsPerRunList = [10755924, 422704, 54287, 4814530, 5425801, 11779739, 36981549, 22712820, 9343259, 5841886, 8387182, 540056, 3832193]

i = 0

for run in runList:

    data = open("fileList"+str(run)+".txt", "r")
    fileList = data.read()

    # print fileList
    print run
    runDirName = "Run"+str(run)
    os.system("mkdir "+runDirName)
    os.system("cat SplitAndRun.py | sed s@THEEVENTS@"+str(eventsPerRunList[i])+"@ >"+runDirName+"/SplitAndRun.py")
    os.system("cp Run.csh "+runDirName)
    os.system("cp runOnCondor "+runDirName)
    # os.system("cat rereco_pp.py | sed s@FILELIST@"+fileList+"@ > "+runDirName+"/rereco_pp.py")
    templateCfg = open("rereco_pp.py").read()
    o = open(runDirName+"/rereco_pp.py","w")
    o.write( re.sub("FILELIST",fileList,templateCfg) )
    o.close()

    os.system("cd "+runDirName+"; python SplitAndRun.py; cd -")

    i += 1
