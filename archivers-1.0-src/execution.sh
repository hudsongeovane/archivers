#!bin/bash
FILE="teste/seq-1to2-2d-2000.txt"
for (( c=1; c<=20; c++ ))
do
./archiver -f $FILE -t 4 -N 100
done

for (( c=1; c<=20; c++ ))
do
./archiver -f $FILE -t 12 -N 100
done

echo "SPEA2 DONE\n"

for (( c=1; c<=20; c++ ))
do
./archiver -f $FILE -t 5 -N 100
done

for (( c=1; c<=20; c++ ))
do
./archiver -f $FILE -t 10 -N 100
done

echo "NSGA2 DONE\n"

for (( c=1; c<=20; c++ ))
do
./archiver -f $FILE -t 6 -N 100
done

for (( c=1; c<=20; c++ ))
do
./archiver -f $FILE -t 13 -N 100
done

echo "AGA DONE\n"

for (( c=1; c<=20; c++ ))
do
./archiver -f $FILE -t 8 -N 100
done

for (( c=1; c<=20; c++ ))
do
./archiver -f $FILE -t 11 -N 100
done

echo "MGA DONE\n"

for (( c=1; c<=20; c++ ))
do
./archiver -f $FILE -t 1 -N 100
done

for (( c=1; c<=20; c++ ))
do
./archiver -f $FILE -t 9 -N 100
done

for (( c=1; c<=20; c++ ))
do
./archiver -f $FILE -t 7 -N 100
done

for (( c=1; c<=20; c++ ))
do
./archiver -f $FILE -t 14 -N 100
done
