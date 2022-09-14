#!/bin/bash
#SBATCH -p cca
#SBATCH -t 10:00:00
#SBATCH -N 1 -n 1
#SBATCH --mem=60G

module purge
module load slurm
module load gcc
module load python3

check_flag() {
	FLAG=$(awk -v n1="$1" -v n2="$2" 'BEGIN {printf ""(n1==n2?"0":"1")"", n1, n2}')
}

#specify dump range and step
DUMPSTART=$1
DUMPSTOP=$2
DUMPSTEP=$3
INC=$4

rm M.txt

M1=1e+20
M2=1e+30

echo $M1 >> M.txt
echo $M2 >> M.txt

MUNIT=$(python3 binsearch.py $DUMPSTART $DUMPSTOP $DUMPSTEP $INC)

NODES=1

#set inclination in model.in
sed -i '/INCLINATION (/s/.*/INCLINATION (deg) '$INC'/' model.in

FLAG=1

while (($FLAG))
do
	echo $MUNIT
	sed -i '/M_UNIT (/s/.*/M_UNIT (g) '$MUNIT'/' model.in	

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
