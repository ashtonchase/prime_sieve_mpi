#! /bin/sh
clear

FILENAME=prime
module load openmpi
g++ ./$FILENAME.cpp -o $FILENAME.o 


    for NUM in 100 1000 10000 100000 1000000 10000000 100000000 1000000000
    do
	
	 $FILENAME.o $NUM | tee $FILENAME.txt
    done


printf "\n\n================================================\n\n"




