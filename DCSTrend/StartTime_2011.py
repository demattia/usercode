#!/usr/bin/python

import time
import os

print "Extracting start and end time from the file"

def makeJsonFile(inputFileName, outputFileName):
    outputFileRun = open(outputFileName, "w")
    outputFileRun.write("{\n")
    outputFileRun.write("    label: 'run',\n")
    outputFileRun.write("    data: [")

    first = True
    file = open(inputFileName)
    start = False
    for line in file:
        if line.find("run") != -1:
            start = True
        elif start and len(line) != 1:
            splittedLine = line.split()
            runNumber = splittedLine[0]

            # results = os.system("./das_client.py --query=\"run = "+str(runNumber)+" | grep run.run_number,run.create_time, run.end_time,run.duration\" --format=plain")
            # print results
            runTimes = os.popen("./das_client.py --query=\"run = "+str(runNumber)+" | grep run.run_number,run.create_time, run.end_time,run.duration\" --format=plain").readlines() 
            print "runTimes =", runTimes

            if len(runTimes) == 0:
                continue

            print "start =", runTimes[0].split(" ")[1].split(".")[0]
            print "end =", runTimes[0].split(" ")[2].split(".")[0]
            startTime = time.mktime(time.strptime(runTimes[0].split(" ")[1].split('.')[0], "%Y-%m-%dT%H:%M:%S"))
            endTime = time.mktime(time.strptime(runTimes[0].split(" ")[2].split('.')[0], "%Y-%m-%dT%H:%M:%S"))
            
            # startTime = splittedLine[1]
            # endTime = splittedLine[2]
            print "run number =", runNumber, "start time =", startTime, "end time =", endTime
            if first:
                first = False
            else:
                outputFileRun.write(", ")

            # print "starttime =", startTime, "(int(startTime)<<16) =", (int(startTime)<<16)
            # print "time =", time.ctime(float(startTime)) #, "and =", time.ctime(float(startTime)<<16)

            outputFileRun.write("["+str(int(startTime)*1000-2000)+", 0, "+str(runNumber)+"], ")
            outputFileRun.write("["+str(int(startTime)*1000)+", 20000, "+str(runNumber)+"], ")
            outputFileRun.write("["+str(int(endTime)*1000)+", 20000, "+str(runNumber)+"], ")
            outputFileRun.write("["+str(int(endTime)*1000+2000)+", 0, "+str(runNumber)+"]")

    outputFileRun.write("]\n}")


makeJsonFile("times_2011.txt", "full_run_2011.js")
# makeJsonFile("times_month.txt", "oneMonth_run.js")

