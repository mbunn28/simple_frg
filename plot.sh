#!/bin/bash

count=15
for (( i = 0; i < $count; ++i )); do
    echo $( echo "-0.5 * $i / $count" | bc -l ) >> tp_vals
    julia -t 2 src/include.jl --nk 16 --nkf 19 --tp $( echo "-0.5 * $i / $count" | bc -l ) --mu $( echo "2 * $i / $count " | bc -l) --U 3 >> lambda_vals
done

# python3 src/plot.py

