#!/bin/sh
#Author domenico.giordano@cern.ch
 
function usage(){
    echo -e "\n[usage] TACAnalysis [options]"
    echo -e " -help  this message"
    echo -e " -fileType=<RU or EDM> (default is ${default_fileType} )"
    echo -e " -tagPN=<tag for PedNoise> (default is ${default_tagPN} )"
    echo -e " -tagCab=<tag for cabling> (default is ${default_tagCab} )"
    echo -e " -InputFilePath=<path>"
    echo -e " -TestArea=<path> (default is $test_area )"
    echo -e " -Flag=<a flag>"
    echo -e " -CondDb=<sqlite>, <devdb10>, <orcon>, <orcoff>, <frontier> (default is ${default_CondDb})"
    echo -e " -sqliteDb=<dbfile> (needed for CondDb=sqlite - default is /tmp/$USER/dummy.db)"
    echo -e " -sqliteCatalog=<dbcatalog> (needed for CondDb=sqlite - default is /tmp/$USER/dummy.db )"
    echo -e " -castor=<file name, or regular-expression> (to get input files from castor)"

    
    echo -e "\nEXAMPLES:"
    echo -e "\n\tSingle Local File access"
    echo -e "\n\t\t./TACAnalysis.sh -CondDb=devdb10 -tagPN=SiStripPedNoise_TOB_v1 -tagCab=SiStripCabling_TOB_v1 -InputFilePath=/storage/TOB/run/RU0002048_000.root -fileType=RU -Flag=Run2048"  

    echo -e "\n\tMultiple Local Files access"
    echo -e "\n\t\t./TACAnalysis.sh -CondDb=devdb10 -tagPN=SiStripPedNoise_TOB_v1 -tagCab=SiStripCabling_TOB_v1 -InputFilePath=/storage/TOB/run/RU0002048_\*.root -fileType=RU -Flag=Run2048"  

    echo -e "\n\tSingle Castor File access"
    echo -e "\n\t\t RU files"
    echo -e "\n\t\t./TACAnalysis.sh -CondDb=devdb10 -tagPN=SiStripPedNoise_TOB_v1 -tagCab=SiStripCabling_TOB_v1 -InputFilePath=/castor/cern.ch/cms/testbeam/TAC/TOB/run -castor=RU0002048_000.root -fileType=RU  -Flag=Run2048"  

    echo -e "\n\t\t EDM files"
    echo -e "\n\t\t./TACAnalysis.sh -CondDb=devdb10 -tagPN=SiStripPedNoise_TOB_v1 -tagCab=SiStripCabling_TOB_v1 -InputFilePath=/castor/cern.ch/cms/store/TAC/TOB/edm_2007_01_20 -castor=EDM0002048_000.root -fileType=EDM  -Flag=Run2048"  

    echo -e "\n\tMultiple Castor Files access (using regular expressions)"

    echo -e "\n\t\t RU files"
    echo -e "\n\t\t./TACAnalysis.sh -CondDb=devdb10 -tagPN=SiStripPedNoise_TOB_v1 -tagCab=SiStripCabling_TOB_v1 -InputFilePath=/castor/cern.ch/cms/testbeam/TAC/TOB/run -castor='RU0002048' -fileType=RU -Flag=Run2048"  

    echo -e "\n\t\t EDM files"
    echo -e "\n\t\t./TACAnalysis.sh -CondDb=devdb10 -tagPN=SiStripPedNoise_TOB_v1 -tagCab=SiStripCabling_TOB_v1 -InputFilePath=/castor/cern.ch/cms/store/TAC/TOB/edm_2007_01_20 -castor=EDM0002048' -fileType=EDM -castor -Flag=Run2048"  
    
    echo
    exit
}


function getLocalRunList(){      
#Create input file list
    
    inputfilenames=""
    for file in `ls ${InputFilePath}`
      do
      [ ! -e $file ] && continue
      inputfilenames="${inputfilenames},\"file:$file\""
    done
    
    inputfilenames=`echo $inputfilenames | sed -e "s@,@@"`
    echo $inputfilenames
}

function getCastorRunList(){      
#Create input file list
    inputfilenames=""
    for file in `nsls ${InputFilePath} | grep -E $castor 2> /dev/null`
      do
      inputfilenames="${inputfilenames},\"castor:${InputFilePath}/$file\""
    done
    
    inputfilenames=`echo $inputfilenames | sed -e "s@,@@"`
    echo $inputfilenames
}

function getRunList(){    

    if [ "$castor" != "0" ]; then
        getCastorRunList
    else
        getLocalRunList 
    fi
}

    function getParameter(){
    what=$1
    shift
    where=$@
    if [ `echo $where | grep -c "\-$what="` = 1 ]; then
        eval $what=`echo $where | awk -F"${what}=" '{print $2}' | awk '{print $1}'`
    elif [ `echo $where | grep -c "\-$what"` = 1 ]; then
	eval $what=1
    else
	let c=$#-1
	shift $c
	eval $what=$1
    fi
}

#################
## MAIN
#################

default_fileType=EDM
default_tagCab=SiStripCabling_TIBD_v1
default_tagPN=SiStripPedNoise_TIBD_v1
default_CondDb=frontier
test_area=/tmp/$USER/TACAnalysis

[ `echo $@ | grep -c "\-help"` = 1 ] && usage;

getParameter InputFilePath $@ .
getParameter TestArea      $@ ${test_area}
getParameter Flag          $@ ""
getParameter CondDb        $@ ${default_CondDb}
getParameter sqliteDb      $@ ${TestArea}/dummy.db
getParameter sqliteCatalog $@ ${TestArea}/dummy.xml
getParameter castor        $@ 0
# Added parameters
# De Mattia 24/1/2007
getParameter tagPN         $@ ${default_tagPN}
getParameter tagCab        $@ ${default_tagCab}
getParameter fileType      $@ ${default_fileType}
# -------------------

[ ! -e ${TestArea} ] && mkdir -p ${TestArea}

frontierFlag=NoFrontier
if [ "$CondDb" == "sqlite" ] && [ "$sqliteDb" != "" ] && [ "$sqliteCatalog" != "" ]; 
    then
    DBfile="sqlite_file:${sqliteDb}"
    DBcatalog="file:${sqliteCatalog}"
elif [ "$CondDb" == "devdb10" ];  then
    DBfile="oracle://devdb10/CMS_COND_STRIP"
    DBcatalog="relationalcatalog_oracle://devdb10/CMS_COND_GENERAL"
elif [ "$CondDb" == "orcon" ]; then
    DBfile="oracle://orcon/CMS_COND_STRIP"
    DBcatalog="relationalcatalog_oracle://orcon/CMS_COND_GENERAL"
elif [ "$CondDb" == "orcoff" ]; then
    DBfile="oracle://cms_orcoff_int2r/CMS_COND_STRIP"
    DBcatalog="relationalcatalog_oracle://cms_orcoff_int2r/CMS_COND_GENERAL"
elif [ "$CondDb" == "frontier" ]; then
    DBfile="frontier://cms_conditions_data/CMS_COND_STRIP"
    DBcatalog=""
    frontierFlag=WithFrontier
else
    echo -e "\nERROR: Wrong options"
    usage
fi

echo -e "\n -InputFilePath=$InputFilePath"
echo -e " -TestArea=$TestArea"
echo -e " -Flag=$Flag"
echo -e " -tagPN=$tagPN"
echo -e " -tagCab=$tagCab"
echo -e " -castor=$castor"
echo -e " -CondDb=$CondDb"
echo -e " -fileType=$fileType"
if [ "$CondDb" = "sqlite" ]; then
    echo -e " -sqliteDb=${sqliteDb}"
    echo -e " -sqliteCatalog=${sqliteCatalog}"
fi
echo " "


cfg_file=${TestArea}/TACAnalysis_${Flag}.cfg
#OUTFILE_CLUSTER=${TestArea}/recocluster_${Flag}.root
OUTFILE_FULL_RECONSTRUCTION=${TestArea}/recofull_${Flag}.root

echo ${cfg_file}

export CORAL_AUTH_PATH=/afs/cern.ch/cms/DB/conddb

eval `scramv1 runtime -sh`

inputfilelist=`getRunList`

[ "$inputfilelist" == "" ] && echo "No file exists for the specified path" && exit

InputSource=PoolSource
[ "$fileType" == "RU" ] && InputSource=TBRUInputSource

templateFile=${CMSSW_BASE}/src/AnalysisExamples/SiStripDetectorPerformance/data/template_TAC_processing.cfg
cat ${templateFile} | sed -e "s#insert_DBfile#$DBfile#" -e "s#insert_DBcatalog#$DBcatalog#" -e "s#insert_input_file_list#$inputfilelist#" -e "s#OUTFILE_CLUSTER#${OUTFILE_CLUSTER}#" -e "s#OUTFILE_FULL_RECONSTRUCTION#${OUTFILE_FULL_RECONSTRUCTION}#" \
-e "s#insert_InputSource#$InputSource#" \
-e "s#insert_tagPN#${tagPN}#g"  -e "s#insert_tagCab#${tagCab}#g" -e "s@#${frontierFlag}@@"> ${cfg_file}

echo "cmsRun ${cfg_file}"

cmsRun ${cfg_file} > ${TestArea}/TACAnalysis_${Flag}.out

echo -e "\nlog file " ${TestArea}/TACAnalysis_${Flag}.out
echo
#echo -e "\nroot files in \n\t\t ${OUTFILE_CLUSTER} \n\t\t ${OUTFILE_FULL_RECONSTRUCTION}\n"
echo -e "\nroot files in \n\t\t ${OUTFILE_FULL_RECONSTRUCTION}\n"
