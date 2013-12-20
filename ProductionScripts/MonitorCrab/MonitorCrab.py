#!/bin/env python

import os, sys
import subprocess

# if len(sys.argv) != 1:
#     print "Please, provide the dir to check"


def checkAndResubmit(crabDir):
    print "Checking status and resubmitting failed jobs for:", crabDir
    # p = subprocess.Popen(["crab", "-status -getoutput -c "+crabDir], stdout=subprocess.PIPE)
    # out, err = p.communicate()
    # print out
    os.system("crab -status -getoutput -c "+crabDir+" >&! monitorTmpLog.txt")

    p = subprocess.Popen(["crab", "-status -c "+crabDir], stdout=subprocess.PIPE)
    out, err = p.communicate()
    # print out

    totalJobs = 0
    jobsWithExitCodeZero = 0
    jobsToResubmit = []
    takeNext = False
    created = False
    createdList = ""
    running = False
    runningList = ""
    cancelled = False
    # cancelledList = ""
    takeNextToNext = 0
    for line in out.split('\n'):
        if line.find("Total") != -1:
            totalJobs = int(line.split()[1])
        # print line
        if takeNextToNext == 1:
            takeNext = True
        if takeNextToNext > 0:
            takeNextToNext -= 1
        if created:
            print line
            createdList += line
            created = False
        if running:
            print line
            runningList += line
            running = False
        if cancelled:
            print line
            os.system("crab -kill "+line.split(':')[1].strip()+" -c "+crabDir)
            # cancelledList += line
            cancelled = False
        if takeNext:
            # print line
            # if jobsToResubmit != "": jobsToResubmit += ','
            # jobsToResubmit += line.split(':')[1].strip()
            jobsToResubmit += [line.split(':')[1].strip()]
            takeNext = False
        if line.find('>>>>') != -1 and line.find('Retrieved') == -1 and line.find('Submitted') == -1 and line.find('Done') == -1:
            if line.find('Created') != -1:
                print "Warning: the following jobs are in status created"
                created = True
            elif line.find('Running') != -1:
                print "The following jobs are in status running"
                running = True
            elif line.find('Cancelled') != -1:
                print "the following jobs are Cancelled and will be resubmitted"
                takeNext = True
                cancelled = True
            elif line.find('Aborted') != -1:
                print "the following jobs are Aborted and will be resubmitted"
                takeNextToNext = 2
            else:
                print line.split(':')
                exitCode = line.split(':')[1].strip()
                if exitCode != "0":
                    print "The following jobs have exitCode =", exitCode, "and will be resubmitted"
                    takeNext = True
                else:
                    jobsWithExitCodeZero = int(line.split()[1])

    if len(jobsToResubmit) == 0:
        print "No job to resubmit at this time"
    else:
        print "\nResubmitting all jobs:", jobsToResubmit

        print jobsToResubmit

    for element in jobsToResubmit:
        print "Resubmitting jobs:", element
        # p = subprocess.Popen(["crab", "-resubmit "+element+" -c "+crabDir], stdout=subprocess.PIPE)
        # out, err = p.communicate()
        # print "out: ", out
        os.system("crab -resubmit "+element+" -c "+crabDir)

    print ""
    print "Summary for "+crabDir+":"
    print "--------\n"
    print "Total number of jobs:", totalJobs
    print "Number of jobs in status done (exit code 0):", jobsWithExitCodeZero
    print "Jobs in status created:"
    print createdList
    print "Jobs in status running:"
    print runningList
    print ""
    if totalJobs == 0:
        print "Error: 0 jobs submitted!"
    else:
        print "Percentange of successfully completed jobs: "+str(int(float(jobsWithExitCodeZero)/float(totalJobs)*100))+"%"
    if totalJobs == jobsWithExitCodeZero:
        print "All jobs completed sucessfully!"
    # Leave an empty line at the end
    print ""


# Running
while True:
    # checkAndResubmit("crabDir/")
