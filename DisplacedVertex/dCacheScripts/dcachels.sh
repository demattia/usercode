#!/bin/sh

subDir=""

if [ "$#" -eq 1 ]; then
    subDir=$1
    # echo "number of parameters = $#"
elif [ "$#" -gt 1 ]; then
    echo "too many parameters. Accepted 1, passed $#"
    exit
fi

echo "srmls -2 \"srm://srm-dcache.rcac.purdue.edu:8443/srm/managerv2?SFN=/store/user/demattia/${subDir}\""
srmls -2 "srm://srm-dcache.rcac.purdue.edu:8443/srm/managerv2?SFN=/store/user/demattia/${subDir}"
