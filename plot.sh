#!/bin/sh

filename=$(basename "$1")
extension="${filename##*.}"
filename="${filename%.*}"

gnuplot -p -e "set terminal pdf;set output '$filename.pdf';unset key;set yrange [0:*];plot '$1'"
