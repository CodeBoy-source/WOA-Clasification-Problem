#!/bin/bash
if [ $# -lt 1 ]
  then
    echo "[ERROR]: You must specify a {seed} and if you want to 0=normal,1=shuffle,2=balance the data "
    exit
fi
# Get Script Directory to later find the bin path
SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]:-$0}"; )" &> /dev/null && pwd 2> /dev/null; )";

if [ -z "$1" ]
  then
    echo "[ERROR]: Couldn't read the seed value for some reason"
    exit
fi

if [ -z "$2" ]
  then
    echo "[ERROR]: Couldn't read the shuffle value for some reason"
    exit
fi

for z in {100..500..50}
do
    for j in 10 30 50 100 500
    do
        for i in $(seq 0.0 .25 0.75)
        do
            echo "[START-$j w $i]: Doing WOAH search in ionosphere.arff"
            $SCRIPT_DIR/../bin/WOAH ionosphere.arff b g 1 $1 $2 $j $i $z 1
            echo "[START-$j w $i]: Doing WOAH search in parkinson.arff"
            $SCRIPT_DIR/../bin/WOAH parkinsons.arff 1 2 1 $1 $2 $j $i $z 1
            echo "[START-$j w $i]: Doing WOAH search in spectf-heart.arff"
            $SCRIPT_DIR/../bin/WOAH spectf-heart.arff 1 2 1 $1 $2 $j $i $z 1
        done;
        for i in $(seq 1.0 .5 4.0)
        do
            echo "[START-$j w $i]: Doing WOAH search in ionosphere.arff"
            $SCRIPT_DIR/../bin/WOAH ionosphere.arff b g 1 $1 $2 $j $i $z 1
            echo "[START-$j w $i]: Doing WOAH search in parkinson.arff"
            $SCRIPT_DIR/../bin/WOAH parkinsons.arff 1 2 1 $1 $2 $j $i $z 1
            echo "[START-$j w $i]: Doing WOAH search in spectf-heart.arff"
            $SCRIPT_DIR/../bin/WOAH spectf-heart.arff 1 2 1 $1 $2 $j $i $z 1
        done;

        for i in 6
        do
            echo "[START-$j w $i]: Doing WOAH search in ionosphere.arff"
            $SCRIPT_DIR/../bin/WOAH ionosphere.arff b g 1 $1 $2 $j $i $z 1
            echo "[START-$j w $i]: Doing WOAH search in parkinson.arff"
            $SCRIPT_DIR/../bin/WOAH parkinsons.arff 1 2 1 $1 $2 $j $i $z 1
            echo "[START-$j w $i]: Doing WOAH search in spectf-heart.arff"
            $SCRIPT_DIR/../bin/WOAH spectf-heart.arff 1 2 1 $1 $2 $j $i $z 1
        done
    done
done
