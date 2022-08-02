#!/bin/bash

echo "RAPTOR setup script"
echo "seting up code for" $1

echo " "

if test -f "raptor_model.c"; then
    echo "raptor_model.c exists. and would be overwritten."
    echo "Only use this script in empty directories."
    echo "Aborting..."
    exit
fi

if [ "$RAPTOR" == "" ] ;
then
	echo "Code dir is empty, set RAPTOR enviroment variable first"
	echo "export RAPTOR=/path/to/code"
	echo "Aborting..."
	exit
else
	echo "Code located at $RAPTOR"
	echo " "
fi

cp $RAPTOR/run/* .

echo "Copying files"
if [ "$1" == "bhac" ] ;
then
	cp -v -n $RAPTOR/model/bhac/* .

	echo "done copying"

elif  [ "$1" == "harm3d" ] ;
then
        cp -v -n $RAPTOR/model/harm3d/* .
else
 echo "no grmhd model picked. Aborting."
 exit
fi

if [ "$1" == "bhac" ] ;
then
	if [ "$2" == "mks" ] ;
	then
		sed -i 's/#define metric (CKS)/#define metric (MKSBHAC)/' parameters.h
	fi
fi

echo " "
echo "Succesfully initialized."
echo " "

echo "Start compilation"
echo " "

make all

if [ $? -ne 0 ]; 
then 
    echo "Compilations unsuccesful" 
    exit
else 
    echo "Compilations succesful" 
fi

echo " "
echo "Deleting setup.sh to ensure you dont overwrite model files from src"

echo " "
echo "Done initialization, use make all to recompile."
