#!/bin/bash

nk=12 #must be even
nkf=25 #must be odd
count=5
for (( i = 0; i < $count; ++i )); do
    echo $( echo "-0.05 * $i / $count - 0.4" | bc -l ) >> tp_vals
    # echo $( echo "-0.4 - 0.2 * $i /$count" | bc -l ) >> mu_vals
    julia -t 2 src/include.jl --nk $nk --nkf $nkf --tp $( echo "-0.05 * $i / $count - 0.4" | bc -l ) --mu $( echo "-0.2 * $i / $count - 1.6" | bc -l ) --U 3 >> lambda_vals
done

# python3 src/plot.py