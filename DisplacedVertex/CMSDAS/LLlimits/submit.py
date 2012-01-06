import sys,os

LeptonType = sys.argv[1]

for mass in range(20,505,5):

	jobfilename = str(mass)+"job"
	job = open(jobfilename,"write")
	job.write("source /afs/cern.ch/user/p/plujan/public/LLlimits/setup.sh\n")
	job.write("python main.py "+str(mass)+ " " + LeptonType +" > /dev/null 2>&1")
	os.system("chmod +x "+jobfilename)
	os.system("bsub -q cmscaf1nh  "+jobfilename)

