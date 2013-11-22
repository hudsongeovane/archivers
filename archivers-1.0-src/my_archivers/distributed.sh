#!bin/bash
FILE="teste/seq-smallPF-2d-10000.txt"
for (( c=1; c<=20; c++ ))
do
./archive $FILE distributed -n 100 -hide -trash
done
