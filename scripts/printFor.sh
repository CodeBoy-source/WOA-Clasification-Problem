#!/bin/bash
for i in {100..500..50}
do
    echo "Results for population of $i"
    ./bin/getMax 5 1 $i
done
