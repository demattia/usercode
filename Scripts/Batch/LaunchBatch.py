#!/usr/bin/python
""" This script uses a cfg file to generate and submit the required jobs
to the batch system. It must be therefore used on lxplus machines.
It performs a check on the existence of the cmssw cfg file.
"""

import sys
import os

print 'cfgFileName =',
if len(sys.argv) != 2:
    print "No cfg name specified, using default batch.cfg"
    print
    f=open('batch.cfg')
else:
    print "Using cfg",
    print sys.argv[1]
    print
    f=open(sys.argv[1])

for line in f:
    # Process only lines that do not start with # and are not empty
    if not(line.startswith("#")) and line != "\n":
        tempTuple = line.split("=")
        varName = tempTuple[0].strip()
        varValue = tempTuple[1].strip()
        # print "varName = " + varName
        if varName == 'cfgName':
            # print "cfg \"" + varValue + "\" found"
            try:
                cfgFile = varValue.split(',')
            except IOError:
                print "ERROR: cfg file \"" + varValue + "\" not found, check",
                print sys.argv[1]
                exit()
        if varName == 'cmsswDir':
            # print "eval scram will be done in " + varValue
            cmsswDir = varValue
        if varName  == 'totalNumberOfEvents':
            totalNumberOfEvents = varValue.split(',')
            # print totalNumberOfEvents
        if varName == 'eventsPerJob':
            eventsPerJob = varValue.split(',')
            # print eventsPerJob
        if varName == 'outFileName':
            outFileName = varValue.split(',')
        if varName == 'outDir':
            outDir = varValue.split(',')
        if varName == 'workingDir':
            workingDir = varValue.split(',')
        if varName == 'randomSeeds':
            randomSeeds = varValue.split(',')
        if varName == 'useSkipEvents':
            useSkipEvents = varValue.split(',')
        if varName == 'pythonCfg':
            pythonCfg = varValue.split(',')

# print "totalNumberOfEvents = " + str(totalNumberOfEvents)
# print "eventsPerJob = " + str(int(eventsPerJob[0]))

# Perform sed-like operations on the cfgFile

typeNum = len(outFileName)

# Check if the number of variables passed is consistent

if (not(len(outDir) == typeNum and len(workingDir) == typeNum and len(totalNumberOfEvents) == typeNum and len(eventsPerJob) == typeNum)):
    print "Number of arguments given in at least one variable not consistent with the others,"
    print "please check that all the multivariables lines in the cfg have the same number of fields separated by commas."
    print "len(outFileName) = " + str(len(outFileName))
    print "len(outDir) = " + str(len(outDir))
    print "len(workingDir) = " + str(len(workingDir))
    print "len(totalNumberOfEvents) = " + str(len(totalNumberOfEvents))
    print "len(eventsPerJob) = " + str(len(eventsPerJob))
    exit()

# Loop on all the types, create subdirs with files to run and submit the jobs
for type in range(typeNum):

    useSkipEventsType = True
    try:
        if ( useSkipEvents[type] == "False" ):
            useSkipEventsType = False
    except (IndexError):
        print "No useSkipEvents specified, using default value True"

    # If all the events are required, create a single job with maxEvents = -1
    if ( totalNumberOfEvents[type].strip() == "-1" ):
        jobsNum = 1
        eventsPerJob[type] = -1
    # Else split the jobs
    else:
        jobsNum = int(totalNumberOfEvents[type])/int(eventsPerJob[type])
    print
    print "Total number of events = " + str(int(totalNumberOfEvents[type]))
    print "Events per job = " + str(eventsPerJob[type])
    print "Creating and submitting " + str(jobsNum) + " jobs"

    skipEvents = 0

    print
    print "Creating job files, any previously existing file will be overwritten"
    print
    os.system("mkdir -p " + workingDir[type])

    for i in range(jobsNum):

        # outFile = str((cfgFile[type].strip()).split('.')[0]).split('/')
        outFile = str((cfgFile[type].strip())).split('.')
        outFile = outFile[len(outFile)-2].split('/')
        outFile = outFile[len(outFile)-1]
        if( not pythonCfg ):
            outFile = open(workingDir[type].strip()+"/"+outFile+"_"+str(i)+".cfg", 'w')
        else:
            outFile = open(workingDir[type].strip()+"/"+outFile+"_"+str(i)+".py", 'w')

        skipEventsFound = False

        # To find maxEvents value, even if on different lines
        maxEventsFound = False

        # print cfgFile[type].strip()
        for s in open(cfgFile[type].strip()):

            if( not(s.startswith("#")) ):
                # Set the maxEvents variable for the job
                if( s.find("maxEvents") != -1 and not pythonCfg ):
                    temp = s.split("=")
                    temp = temp[2].split("}")[0].strip()
                    # print "temp = " + temp
                    outFile.write(s.replace( temp, str(eventsPerJob[type]) ))

                # Set the maxEvents variable for the job
                if( (s.find("maxEvents") != -1 or not maxEventsFound) and pythonCfg ):
                    if( s.find("input") != -1 ):
                        temp = s.split("(")
                        temp = temp[2].split(")")[0].strip()
                        maxEventsFound = True
                        # print "temp = " + temp
                        outFile.write(s.replace( temp, str(eventsPerJob[type]) ))
                    else:
                        outFile.write(s);
                # Set the skipEvents variable
                elif( s.find("skipEvents") != -1 and useSkipEventsType == True ):
                    temp = s.split("=")
                    temp = temp[1].strip()
                    # print temp
                    outFile.write(s.replace( temp, str(skipEvents) ))
                    skipEventsFound = True
                # Change the output file name
                elif( s.find(outFileName[type].strip()) != -1 ):
                    temp = (outFileName[type].split('.')[0]).strip()
                    outFile.write(s.replace( temp, temp + "_" + str(i) ) )
                else:
                    # Loop on all the random seeds and change them
                    seedFound = False
                    for seed in randomSeeds:
                        #if (s.find(seed) != -1 and (s.find(seed+str('.')) == -1 ) and (s.find("replace") == -1 )):
                        if (s.find(seed) != -1):
                            temp = s.split("=")
                            try:
                                if( not pythonCfg ):
                                    seedValue = int(temp[1].strip())
                                if( pythonCfg ):
                                    seedValue = int(temp[1].split("(")[1].split(")")[0].strip())
                            # If the name of the seed matches another name there could be two cases:
                            # 1) There is no value beyond the "=" (or no equal at all) -> IndexError is thrown
                            # 2) The value past the "=" is not convertible to an integer -> ValueError is thrown
                            # In both cases the unchanged line is written.
                            # ATTENTION: this script cannot distinguish between real seeds and no seeds. If the
                            # expression seedName = seedValue is matched it will change the second number.
                            except (ValueError, IndexError):
                                break
                            # Otherwise the seed value is changed
                            else:
                                outFile.write(s.replace( str(seedValue), str(seedValue + i) ))
                                seedFound = True
                    if (not(seedFound)):
                        outFile.write(s)
            else:
                outFile.write(s)

        if skipEventsFound == False and useSkipEvents == True:
            print 'No skipEvents field found in the cfg file, please add it to the source module'
            exit()

        batchFileName = workingDir[type].strip() + "/" + "batch_" + str(i) + ".csh"
        batchFile = open(batchFileName, 'w')
        for s in open('batch_test.csh'):
            # Set the cmssw dir where to do eval scram
            if( s.find("cmsswDir") != -1 ):
                batchFile.write(s.replace('cmsswDir', cmsswDir.strip()))
            elif( s.find("cfgFile") != -1 ):
                tempCfgFileName = cfgFile[type].split('/')
                tempCfgFileName = (tempCfgFileName[len(tempCfgFileName)-1].split('.')[0]).strip()
                currentDir = os.getcwd()
                if( not pythonCfg ):
                    batchFile.write(s.replace('cfgFile', currentDir.strip() + "/" + workingDir[type].strip() + "/"+ tempCfgFileName+'_'+str(i)+'.cfg'))
                else:
                    batchFile.write(s.replace('cfgFile', currentDir.strip() + "/" + workingDir[type].strip() + "/"+ tempCfgFileName+'_'+str(i)+'.py'))
                # Get the current working dir
            elif( s.find("outFileName") != -1 ):
                temp = outFileName[type].split('.')
                batchFile.write(s.replace('outFileName', temp[0].strip()+"_"+str(i)+".root").replace('outDir', outDir[type].strip()))
            else:
                batchFile.write(s)

        # Close the file here or will not be executable later (busy)
        batchFile.close()

        skipEvents += int(eventsPerJob[type])

        # Make the batch script executable and run it
        os.system("chmod 777 " + batchFileName)
        #os.system("bsub -R \"pool>40\" -q 8nh -J " + tempCfgFileName + "<" + batchFileName)

    # end of loop on jobs

# end of loops on types

print
print "The batch jobs will fill your home dir with LSF* directories containing stdout and stderr files."
print "You can use this command \"find . -name \"LSF*\" -exec rm -rf {} \;\" to remove all of them"
print

