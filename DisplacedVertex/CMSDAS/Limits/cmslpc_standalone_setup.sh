#!/bin/sh
#
# This script will set up the specified ROOT release
# and configure the corresponding environment on lxplus.
# Version for bash-like shell.
#
# 2011 Gena Kukartsev
#
# Usage:
#        source cmslpc_standalone_setup.sh
#


echo ''
echo 'Setting up python, ROOT and PyROOT'
echo ''

#arch=slc5_ia32_gcc434
arch=slc5_amd64_gcc434

export CMS_PATH=/uscmst1/prod/sw/cms

#export PYTHONDIR=/uscmst1/prod/sw/cms/slc5_ia32_gcc434/external/python/2.6.4-cms8
export PYTHONDIR=/uscmst1/prod/sw/cms/slc5_amd64_gcc434/external/python/2.6.4-cms15

#export PATH=${PYTHONDIR}/bin:/uscms/home/kukarzev/nobackup/root/root_v5.28.00b/bin:$PATH
#export PATH=${PYTHONDIR}/bin:/uscms/home/kukarzev/nobackup/root/root_v5.30.00/bin:$PATH
#export PATH=${PYTHONDIR}/bin:/uscms/home/kukarzev/nobackup/root/root_v5.30.02_amd64/bin:$PATH
export PATH=${PYTHONDIR}/bin:/uscms/home/kukarzev/nobackup/root/root_v5.32.00/bin:$PATH

#export ROOTSYS=/uscms/home/kukarzev/nobackup/root/root_v5.28.00b
#export ROOTSYS=/uscms/home/kukarzev/nobackup/root/root_v5.30.00
#export ROOTSYS=/uscms/home/kukarzev/nobackup/root/root_v5.30.02
#export ROOTSYS=/uscms/home/kukarzev/nobackup/root/root_v5.30.02_amd64
export ROOTSYS=/uscms/home/kukarzev/nobackup/root/root_v5.32.00

export PYTHONPATH=${ROOTSYS}/lib:${PYTHONPATH}

export LD_LIBRARY_PATH=${PYTHONDIR}/lib:${CMS_PATH}/$arch/external/gcc/4.3.4/lib64:${CMS_PATH}/$arch/external/gcc/4.3.4/lib:${ROOTSYS}:${ROOTSYS}/lib:${LD_LIBRARY_PATH}

export ROOT_INCLUDE=${ROOTSYS}/include
