#!bin/bash
FILE="teste/seq-1to2-2d-2000.txt"
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
