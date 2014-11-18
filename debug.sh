#!/bin/sh

filename=$(basename "$1")
extension="${filename##*.}"
filename="${filename%.*}"

# gnuplot -p -e "set terminal pdf;set output '$filename.pdf';unset key;set yrange [0:*];plot '$1'"
gnuplot -p -e "set terminal pdf;set output '$filename.pdf';set yrange [0:*]; plot '$1' i 0 u 1:2 w lines title columnheader(1), '$1' i 1 u 1:2 w lines title columnheader(1)"
#gnuplot -p -e "set terminal pdf;set output '$filename.pdf';set yrange [1:*];set logscale y; plot '$1' i 0 u 1:2 w dots title columnheader(1), '$1' i 1 u 1:2 w dots title columnheader(1)"
#gnuplot -p -e "set terminal pdf;set output '$filename.pdf';set yrange [1:*];set logscale y; plot '$1' i 0 u 1:2 w dots title columnheader(1), '$1' i 1 u 1:2 w dots title columnheader(1)"
#gnuplot -p -e "set terminal pdf;set output '$filename.pdf';set yrange [1:1.21e+24];set logscale y; set ytics(1.02e+03,1.05e+06,1.07e+09,1.1e+12,1.13e+15,1.15e+18,1.18e+21,1.21e+24);plot '$1' i 0 u 1:2 title columnheader(1)"
#gnuplot -p -e "set terminal pdf;set output '$filename.pdf';set yrange [1:*];set logscale y; plot '$1' i 0 u 1:2 w dots title columnheader(1), '$1' i 1 u 1:2 w dots title columnheader(1), '$1' i 2 u 1:2 w dots title columnheader(1)"
