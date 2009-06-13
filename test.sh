#!/bin/bash
# Usage: test -b

### Start main procedure

for n in 500 1000 2000 4000 8000
do
	N=`expr $n \* 5`
	if [ $1 = "-b" ]
	then
		./build.sh -n $n -N $N
	else
		echo && echo && echo && echo "Testing with $n initial particles per node" >> test.log
		echo "Testing with $n initial particles per node" >> test.log
		date +"%Y-%m-%d %H:%M" >> test.log
		time mpprun ./gaslaw_${n}_$N >> test.log
	fi
done