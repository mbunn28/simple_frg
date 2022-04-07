#!/bin/bash

count=10
for (( i = 0; i < $count; ++i )); do
    echo $( echo "-0.5 * $i / $count" | bc -l ) >> tp_vals
    julia -t 2 src/include.jl --nk 12 --nkf 15 --tp $( echo "-0.5 * $i / $count" | bc -l ) --mu $( echo "-2 * $i / $count " | bc -l) --U 3 >> lambda_vals
done

# python3 src/plot.py

