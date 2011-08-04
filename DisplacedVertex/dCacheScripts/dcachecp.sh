#!/bin/sh

type=""

if [ "$#" -eq 2 ]; then
    type=$1
    destinationDir=$2
else
    echo "Wrong number of parameters. Accepted 1, passed $#"
    exit
fi

if [[ ${type} =~ '^/' ]]; then
    echo "full path"
else
    echo "local dir"
    fullPath=`pwd`
    type=`echo ${fullPath}/${type}`
fi

for file in ${type}; do
    echo "${file}"
    filename=$(basename ${file})
    echo "lcg-cp --verbose -b -D srmv2 \"file:${file}\" \"srm://srm-dcache.rcac.purdue.edu:8443/srm/managerv2?SFN=/store/user/demattia/${destinationDir}/${filename}\""
    lcg-cp --verbose -b -D srmv2 "file:${file}" "srm://srm-dcache.rcac.purdue.edu:8443/srm/managerv2?SFN=/store/user/demattia/${destinationDir}/${filename}"
done

