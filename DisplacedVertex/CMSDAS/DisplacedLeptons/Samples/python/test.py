import subprocess

# output = subprocess.Popen(['ls', '-l'], stdout=subprocess.PIPE).communicate()[0]
result = subprocess.Popen(['srmls', "-2 -offset 0 \"srm://srm-dcache.rcac.purdue.edu:8443/srm/managerv2?SFN=/store/user/demattia/longlived/CMSSW_4_2_7/QCDmu80/pat/\""], stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()[0]

# print result
for line in result.splitlines():
    print "line = "+line
