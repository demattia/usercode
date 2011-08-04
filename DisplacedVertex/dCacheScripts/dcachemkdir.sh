#!/bin/sh

subDir=""

if [ "$#" -eq 1 ]; then
    subDir=$1
else
    echo "Wrong number of parameters. Accepted 1, passed $#"
    exit
fi

echo "srmmkdir -2 \"srm://srm-dcache.rcac.purdue.edu:8443/srm/managerv2?SFN=/store/user/demattia/${subDir}\""
srmmkdir -2 "srm://srm-dcache.rcac.purdue.edu:8443/srm/managerv2?SFN=/store/user/demattia/${subDir}"

