#!/bin/bash

echo "RAPTOR setup script"
echo "seting up code for" $1

echo " "


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

if [ "$1" == "bhac" ] ;
then
	cp $RAPTOR/model/bhac/* .

elif  [ "$1" == "harm3d" ] ;
then
 cp $RAPTOR/model/harm3d/* .
else
 echo "no grmhd model picked. Aborting."
 exit
fi

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

rm setup.sh

echo " "
echo "Done initialization, use make all to recompile."
