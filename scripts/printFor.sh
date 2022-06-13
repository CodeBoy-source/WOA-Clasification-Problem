#!/bin/bash
for i in {100..500..50}
do
    echo "Results for population of $i"
    ./bin/getMax 5 1 $i
done
for i in 15000
do
    echo "Results for population of $i"
    ./bin/getMax 5 1 $i
done
