#!bin/bash
FILE="teste/seq10000points_1000pf_4D.txt"
for (( c=1; c<=20; c++ ))
do
./archive $FILE ara -n 100 -hide
done


for (( c=1; c<=20; c++ ))
do
./archive $FILE ara -n 100 -hide -trash
done
