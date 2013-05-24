import os

def submit(dirName):
    print "cd "+dirName+"; rm -f condor_rereco*; rm -f selection_test.root; condor_submit condor; cd .."
    os.system("cd "+dirName+"; rm -f condor_rereco*; rm -f selection_test.root; condor_submit condor; cd ..")

dirList = ["Run2012A", "Run2012ARecover", "Run2012B", "Run2012C1", "Run2012DRereco", "Run2012C2", "Run2012D", "BsMC"]
# dirList = ["Run2012CEcalRecover", "Run2012C2", "Run2012D", "BsMC"]
# dirList = ["Run2012A", "Run2012ARecover", "Run2012B", "Run2012C1"]

for name in dirList:
    submit(name)
