#!/bin/sh

mainDir="/uscms_data/d3/demattia/DisplacedLeptons/Trigger/OpenHLT/"
declare -a names=('MH120MFF50' 'MH200MFF20' 'MH200MFF50' 'MH400MFF05' 'MH400MFF20' 'MH400MFF50' 'MH400MFF150' 'MH1000MFF20' 'MH1000MFF50' 'MH1000MFF150' 'MH1000MFF350' 'part2')

num=${#names[@]}
echo "number of samples = ${num}"

i=0

localDir=`pwd`

while [ "$i" -lt "$num" ]; do
    cd ${localDir}
    cp checkOpenHLT.* ${mainDir}${names[$i]}
    cp drawCanvases.C ${mainDir}${names[$i]}
    cp FillPage.sh ${mainDir}${names[$i]}
    cp template_page.html ${mainDir}${names[$i]}
    cd ${mainDir}${names[$i]}
    root -l -b -q checkOpenHLT.C+
    ./FillPage.sh
    rm -rf ${localDir}/${names[$i]}
    mv WebPage ${localDir}/${names[$i]}
    
    let "i+=1"
done
