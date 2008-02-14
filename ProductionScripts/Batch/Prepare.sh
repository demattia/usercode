#!/bin/sh

if [ "$1" == "" ]; then
	echo Provide a bin value, eg 230_300, underscore required
else

	name2=`echo $1 | awk -F_ '{print $1 "-" $2}'`

	echo name2 = $name2

	mkdir QCD_$name2

	cat QCD_170-230/Fast_qcd_170_230_puhl.cfg | sed s/170_230/$1/g | sed s/170-230/$name2/g > QCD_$name2/Fast_qcd_$1_puhl.cfg

	cat QCD_170-230/batch_qcd_170_230_puhl.csh | sed s/170_230/$1/g | sed s/170-230/$name2/g > QCD_$name2/batch_qcd_$1_puhl.csh

	cat QCD_170-230/Generator.sh | sed s/170_230/$1/g | sed s/170-230/$name2/g > QCD_$name2/Generator.sh

fi
