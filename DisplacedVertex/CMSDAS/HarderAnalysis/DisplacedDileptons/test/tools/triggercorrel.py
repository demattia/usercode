#!/usr/bin/python

# small tool to investigate differences between trigger selections
# in Paul Lujan's ASCII files

import os,sys

ntrig=0
trig=[]
trigsums={}
for filename in os.listdir("."):
    if filename.find("triggerResults_")!=0: continue
    print "processing",filename,"..."
    results=open(filename,"r").readlines()
    for line in results:
        if line.find("HLT")>=0: continue
        if len(line.split())==1 and ntrig==0:
            ntrig=int(line.strip())
            if ntrig>10:
                print "error: cannot deal with more than 10 triggers"
                sys.exit(1)
            for i in range(ntrig):
                trig.append(0)
                for k in range(i+1,ntrig):
                    trigsums[i*1000    +k*10  ]=0
                    trigsums[i*1000+100+k*10+1]=0
                    trigsums[i*1000    +k*10+1]=0
                    trigsums[i*1000+100+k*10  ]=0
        if len(line.split())<6: continue
        for i in range(ntrig):
            trig[i]=int(line.split()[3+i])
        for i in range(ntrig):
            for k in range(i+1,ntrig):
                trigsums[i*1000+trig[i]*100+k*10+trig[k]]+=1

print "RESULTS:"
for i in range(ntrig):
    for k in range(i+1,ntrig):
        print "triggers",i,"and",k
        print "  both fired:    %7i"%trigsums[i*1000+100+k*10+1]
        print "  neither fired: %7i"%trigsums[i*1000+k*10]
        print "  only",i,"fired : %7i"%trigsums[i*1000+100+k*10]
        print "  only",k,"fired : %7i"%trigsums[i*1000+k*10+1]
