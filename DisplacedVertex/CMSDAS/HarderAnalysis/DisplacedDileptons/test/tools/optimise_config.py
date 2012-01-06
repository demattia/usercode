#!/bin/env python

import os,sys,time
today=time.strftime("%d-%b-%Y")


def run_with_skipped_line(origcode,line_to_be_skipped,logfile):
    if (line_to_be_skipped>=0):
        message="OPTIMISER: skipping line"+str(line_to_be_skipped+1)+":"+\
                 origcode[line_to_be_skipped].strip()
    else:
        message="OPTIMISER: running full file"

    # create modified config file
    outfile=open("optimiser_test.py","w")
    for i in range(len(origcode)):
        if i!=line_to_be_skipped:
            outfile.write(origcode[i])
    outfile.close()

    # run the modified config file
    os.system("rm -f optimiser_test.log")
    os.popen("cmsRun optimiser_test.py &> optimiser_test.log")

    # interpret output
    resultfile=open("optimiser_test.log","r")
    failure=0
    success=0
    for line in resultfile.readlines():
        if line.find("exception")>=0: failure=1
        if line.find("TrigReport")==0:
            words=line.split()
            if len(words)==7:
                if (words[6]=="out") and (int(words[3])>0) and (words[4]=="0") and (words[5]=="0"):
                    success=1

    # if the job did not finish or caused exceptions, then we are done
    if failure or not success:
        message+=" ==> FAILS"
        return message

    # the PAT creation seems to have worked; but is the analysis still getting the same results?
    os.popen("cmsRun opt_cfg.py &> opt.log")
    if line_to_be_skipped>=0:
        os.system("grep -v "+today+" opt.log > opt2.log")
        diffs=os.popen("diff opt2.log opt_ref2.log").readlines()
        if len(diffs)>0:
            message+=" => DIFFERENT RESULT"
            return message

    message+=" ==> WORKS"
    return message


#***********************************************
#                 MAIN CODE
#***********************************************

# read original config file
origfile=open("HarderAnalysis/DisplacedDileptons/test/makePATtupleReReco_cfg.py","r")
origcode=origfile.readlines()
origfile.close()

# prepare e and mu config files for PAT analysis
os.system("cp HarderAnalysis/DisplacedDileptons/test/main_cfg.py opt_cfg.py")
cfgfile=open("opt_cfg.py","a")
cfgfile.write("\nprocess.source.fileNames=cms.untracked.vstring(\"file:PATtupleReReco_new.root\")\n")
cfgfile.close()

# run the complete file for reference
logfile=open("optimiser.summary","w")
result=run_with_skipped_line(origcode,-1,logfile)
print result
logfile.write(result+"\n")
logfile.close()

# create all the reference files for later
os.system("mv histograms.root histograms_optref.root")
os.system("mv opt.log opt_ref.log")
os.system("grep -v "+today+" opt_ref.log > opt_ref2.log")

# now loop over all lines and remove each one at a time
line_to_be_skipped=-1
while line_to_be_skipped<len(origcode)-1:
    line_to_be_skipped+=1
    line=origcode[line_to_be_skipped].strip()
    if line.find("keep")<0: continue          # special case: optimise only output collections
    if len(line)<5: continue
    if line[0]=="#": continue
    logfile=open("optimiser.summary","a")
    result=run_with_skipped_line(origcode,line_to_be_skipped,logfile)
    print result
    logfile.write(result+"\n")
    logfile.close()

# remove temporary files
os.system("rm optimiser_test.py optimiser_test.log")
os.system("rm opt.log opt2.log opt_ref.log histograms_optref.root opt_cfg.py")
