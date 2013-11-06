#!bin/bash
FILE="teste/seq-smallPF-3d-10000.txt"
for (( c=1; c<=20; c++ ))
do
./archiver -f $FILE -t 5 -N 100
done

for (( c=1; c<=20; c++ ))
do
./archiver -f $FILE -t 10 -N 100
done
