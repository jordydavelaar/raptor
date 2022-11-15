#!/bin/bash
#SBATCH -p cca
#SBATCH -t 10:00:00
#SBATCH -N 1 -n 1
#SBATCH --mem=60G

module load hdf5 gsl
module load python3

check_flag() {
	FLAG=$(awk -v n1="$1" -v n2="$2" 'BEGIN {printf ""(n1==n2?"0":"1")"", n1, n2}')
}

#specify dump range and step
DUMPSTART=$1
DUMPSTOP=$2
DUMPSTEP=$3

NODES=1

maindir=$PWD

#array for variables we want to run
inc=(20)
Rh=(1)

#range for Munit
M1=1e+20
M2=1e+30


function binsearch {


   mkdir $maindir/i$1
   mkdir $maindir/i$1/R$2
   cd $maindir/i$1/R$2

   INC=$1
   RHIGH=$2

   $RAPTOR/setup.sh -c=bhac -m=mks -r=unpol -s=sfc

   cp $maindir/rap.sh $maindir/i$1/R$2/.
   cp $maindir/binsearch.py $maindir/i$1/R$2/. 
   cp $maindir/grid_mks.in $maindir/i$1/R$2/.

   rm M.txt

   echo $M1 >> M.txt
   echo $M2 >> M.txt

   MUNIT=$(python3 binsearch.py $DUMPSTART $DUMPSTOP $DUMPSTEP $INC)

   #set inclination in model.in
   sed -i '/INCLINATION\t(/s/.*/INCLINATION\t(deg)\t'$INC'/' model.in
   sed -i '/R_HIGH\t\t(/s/.*/R_HIGH\t\t(-)\t'$RHIGH'/' model.in

   FLAG=1

   while (($FLAG))
   do
	echo $MUNIT
	sed -i '/M_UNIT\t\t(/s/.*/M_UNIT\t\t(g)\t'$MUNIT'/' model.in	

	# the loop below works, but it seems the waiting does not ... weird
	#for i in $(seq $DUMPSTART $DUMPSTEP $(expr $DUMPSTOP - $DUMPSTEP))
	#do
	#sbatch --wait rap_mult.sh $i $(expr $i + $DUMPSTEP) $DUMPSTEP 90 $MUNIT & 
	#done

	DELTAI=$(expr $DUMPSTOP - $DUMPSTART)

	let DELTAI=$DELTAI/$NODES

	echo $DELTAI

	for((i=0;i<$NODES ;i++))
	do
        	START=$(expr $i \* $DELTAI + $DUMPSTART)
	        STOP=$(expr $START + $DELTAI)
	        echo $START $STOP $DUMPSTEP
	        sbatch --wait rap.sh $START $STOP $DUMPSTEP &
	done

	if(($DUMPSTOP>$STOP))
	then
        	LEFT=$(expr $DUMPSTOP - $STOP )
	        START=$(expr $STOP )
	        STOP_LEFT=$(expr $START + $LEFT )
	        echo $START $STOP_LEFT $DUMPSTEP
	        sbatch --wait rap.sh $START $STOP_LEFT $DUMPSTEP &

	fi

	wait

	MUNIT=$(python3 binsearch.py $DUMPSTART $DUMPSTOP $DUMPSTEP $INC)

	second_to_last=$(tail -n 2 M.txt | head -n 1)
	last=$(tail -n 1 M.txt)

	check_flag $second_to_last $last
	echo "Will continue if not zero: " $FLAG
   done
}

for i in ${inc[@]}; do
  for r in  ${Rh[@]}; do
    binsearch $i $r
  done
done

