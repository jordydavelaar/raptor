#!/bin/bash

echo "RAPTOR setup script"


METRIC="cks"
INT="rk4"
CODE="bhac"
GRID="smr"
RAD="pol"

for arg in "$@"
do
    case $arg in
        -c=*|--code=*)
        CODE="${arg#*=}"
        shift # Remove --cache= from processing
        ;;
	-m=*|--metric=*)
	METRIC="${arg#*=}"
	shift # Remove --cache= from processing
	;;
	-i=*|--int=*)
	INT="${arg#*=}"
        shift # Remove --cache= from processing
        ;;
	-g=*|--grid=*)
	GRID="${arg#*=}"
        shift # Remove --cache= from processing
        ;;
	-r=*|--rad=*)
        RAD="${arg#*=}"
        shift # Remove --cache= from processing
        ;;	

    esac
done

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
if [ "$CODE" == "bhac" ] ;
then
	cp -v -u -n $RAPTOR/model/bhac/* .

	echo "done copying"

elif  [ "$CODE" == "harm3d" ] ;
then
        cp -v -u -n $RAPTOR/model/harm3d/* .
else
 echo "no grmhd model picked. Aborting."
 exit
fi

if [ "$CODE" == "bhac" ] ;
then
	if [ "$METRIC" == "mks" ] ;
	then
		sed -i '/#define metric (/s/.*/#define metric (MKSBHAC)/' parameters.h
	fi
    	if [ "$METRIC" == "cks" ] ;
    	then
	    	sed -i '/#define metric (/s/.*/#define metric (CKS)/' parameters.h
        fi	
fi

if [ "$CODE" == "harm3d" ] ;
then
    	if [ "$METRIC" == "mks" ] ;
    	then
	    	sed -i '/#define metric (/s/.*/#define metric (MKSHARM)/' parameters.h
        fi
fi

if [ "$INT" == "ver" ] ;
then
      	sed -i '/#define int_method (/s/.*/#define int_method (VER)/' parameters.h
fi


if [ "$INT" == "rk2" ] ;
then
      	sed -i '/#define int_method (/s/.*/#define int_method (RK2)/' parameters.h
fi

if [ "$INT" == "rk4" ] ;
then
      	sed -i '/#define int_method (/s/.*/#define int_method (RK4)/' parameters.h
fi

if [ "$INT" == "rk45" ] ;
then
      	sed -i '/#define int_method (/s/.*/#define int_method (RK45)/' parameters.h
fi

if [ "$GRID" == "SMR" ] ;
then

      	sed -i '/#define AMR /s/.*/#define AMR 0/' parameters.h
      	sed -i '/#define SMR /s/.*/#define SMR 1/' parameters.h
fi

if [ "$GRID" == "AMR" ] ;
then
    	sed -i '/#define AMR /s/.*/#define AMR 1/' parameters.h
        sed -i '/#define SMR /s/.*/#define SMR 0/' parameters.h
fi

if [ "$RAD" == "pol" ] ;
then
      	sed -i '/#define POL (/s/.*/#define POL (1)/' parameters.h
fi

if [ "$RAD" == "unpol" ] ;
then
    	sed -i '/#define POL (/s/.*/#define POL (0)/' parameters.h
fi

echo " "
echo "Succesfully initialized."
echo " "

echo "Start compilation"
echo " "

make clean && make all

if [ $? -ne 0 ]; 
then 
    echo "Compilations unsuccesful" 
    exit
else 
    echo "Compilations succesful" 
fi

echo " "
echo "Done initialization, use make all to recompile."

echo "Setup used code $CODE (-c/--code) metric $METRIC (-m/--metric) int $INT (-i/--int) grid $GRID (-g/--grid) rad $RAD (-r/--rad)"
