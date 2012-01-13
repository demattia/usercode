#!/bin/tcsh
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
set arch=slc5_amd64_gcc434

setenv CMS_PATH /uscmst1/prod/sw/cms

#setenv PYTHONDIR /uscmst1/prod/sw/cms/slc5_ia32_gcc434/external/python/2.6.4-cms8
setenv PYTHONDIR /uscmst1/prod/sw/cms/slc5_amd64_gcc434/external/python/2.6.4-cms15

#setenv PATH ${PYTHONDIR}/bin:/uscms/home/kukarzev/nobackup/root/root_v5.28.00b/bin:$PATH
#setenv PATH ${PYTHONDIR}/bin:/uscms/home/kukarzev/nobackup/root/root_v5.30.00/bin:$PATH
#setenv PATH ${PYTHONDIR}/bin:/uscms/home/kukarzev/nobackup/root/root_v5.30.02_amd64/bin:$PATH
setenv PATH ${PYTHONDIR}/bin:/uscms/home/kukarzev/nobackup/root/root_v5.32.00/bin:$PATH

#setenv ROOTSYS /uscms/home/kukarzev/nobackup/root/root_v5.28.00b
#setenv ROOTSYS /uscms/home/kukarzev/nobackup/root/root_v5.30.00
#setenv ROOTSYS /uscms/home/kukarzev/nobackup/root/root_v5.30.02
#setenv ROOTSYS /uscms/home/kukarzev/nobackup/root/root_v5.30.02_amd64
setenv ROOTSYS /uscms/home/kukarzev/nobackup/root/root_v5.32.00

setenv PYTHONPATH ${ROOTSYS}/lib:${PYTHONPATH}

setenv LD_LIBRARY_PATH ${PYTHONDIR}/lib:${CMS_PATH}/$arch/external/gcc/4.3.4/lib64:${CMS_PATH}/$arch/external/gcc/4.3.4/lib:${ROOTSYS}:${ROOTSYS}/lib:${LD_LIBRARY_PATH}

setenv  ROOT_INCLUDE ${ROOTSYS}/include
