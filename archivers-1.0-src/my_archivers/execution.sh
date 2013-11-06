#!bin/bash
FILE="teste/seq2000points_500pf_3D.txt"
for (( c=1; c<=20; c++ ))
do
./archive $FILE ideal -n 100 -hide
done


for (( c=1; c<=20; c++ ))
do
./archive $FILE ideal -n 100 -hide -trash
done

for (( c=1; c<=20; c++ ))
do
./archive $FILE distributed -n 100 -hide
done


for (( c=1; c<=20; c++ ))
do
./archive $FILE distributed -n 100 -hide -trash
done

for (( c=1; c<=20; c++ ))
do
./archive $FILE distance -n 100 -hide
done


for (( c=1; c<=20; c++ ))
do
./archive $FILE distance -n 100 -hide -trash
done

for (( c=1; c<=20; c++ ))
do
./archive $FILE ara -n 100 -hide
done


for (( c=1; c<=20; c++ ))
do
./archive $FILE ara -n 100 -hide -trash
done
