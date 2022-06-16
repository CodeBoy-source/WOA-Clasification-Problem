#!/bin/bash
if [ $# -lt 2 ]
  then
    echo "[ERROR]: You must specify a {seed}, if you want to 0=normal,1=shuffle,2=balance the data and the population size. "
    echo "[EXAMPLE]: ./WOAH.sh 150421 0 10 "
    echo "[DEFAULT VALUES ONLY] For testing lambda and/or delta run manually as for instance: "
    echo "./bin/WOAH parkinsons.arff 1 2 1 150421 0 10 (lambda) (delta)"
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

if [ -z "$3" ]
  then
    echo "[ERROR]: Couldn't read the population size for some reason"
    exit
fi

echo "[START-1]: Doing WOAH search in ionosphere.arff"
$SCRIPT_DIR/../bin/WOAH ionosphere.arff b g 1 $1 $2 $3
echo "[START-2]: Doing WOA search in parkinson.arff"
$SCRIPT_DIR/../bin/WOAH parkinsons.arff 1 2 1 $1 $2 $3
echo "[START-3]: Doing WOA search in spectf-heart.arff"
$SCRIPT_DIR/../bin/WOAH spectf-heart.arff 1 2 1 $1 $2 $3


