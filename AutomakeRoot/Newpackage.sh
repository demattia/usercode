#!/bin/sh

if [ "$1" == "" ]; then
    echo Please, provide a name for the package
    exit
fi

echo Generating package files
cat template_Makefile.am | sed s/hello/$1/g > Makefile.am
cat template_configure.in | sed s/hello/$1/g > configure.in
mkdir -p src
cat template_src_Makefile.am | sed s/hello/$1/g > src/Makefile.am
cat template_hello.cpp | sed s/hello/$1/g > src/$1.cpp

# echo deleting templates
# rm template*

echo executing acloacl
echo
aclocal
echo
echo running autoconf
echo
autoconf
echo
echo files README, AUTHORS, NEWS and ChangeLog created
echo
touch README AUTHORS NEWS ChangeLog
echo
echo creating config.h with autoheader
echo
autoheader
echo
echo running automake ...
echo
automake -a 
echo
echo executing configure ...
echo
./configure
echo
echo making the package ...
echo
make
