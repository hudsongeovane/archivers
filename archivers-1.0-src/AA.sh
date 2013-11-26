#!bin/bash
FILE="teste/seq10000points_1000pf_4D.txt"
for (( c=1; c<=20; c++ ))
do
./archiver -f $FILE -t 7 -N 100
done

for (( c=1; c<=20; c++ ))
do
./archiver -f $FILE -t 14 -N 100
done
