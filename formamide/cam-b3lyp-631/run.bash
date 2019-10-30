#!/bin/bash

g09 td1.com
g09 td2.com
g09 td3.com
g09 td4.com
g09 td5.com
g09 td6.com
g09 td7.com
g09 td8.com

g09 td1-opt.com
g09 td2-opt.com
g09 td3-opt.com
g09 td4-opt.com
g09 td5-opt.com
g09 td6-opt.com
g09 td7-opt.com
g09 td8-opt.com

formchk freq.chk
for i in {1..8}; do 
    formchk td${i}.chk
    cp td${i}.fchk td.fchk
    cp td${i}-opt.fchk td-opt.fchk
    ./calc_deltas.bin 6 1
    cp fort.77 disp${i}.dat
done

for i in {1..8}; do 
    grep -A12 LMXg disp${i}.dat > tmpopt${i}
    grep -A12 gradient disp${i}.dat > tmpgrad${i}
done

paste tmpopt1 tmpopt2 tmpopt3 tmpopt4 tmpopt5 tmpopt6 tmpopt7 tmpopt8 > dispopt.dat
paste tmpgrad1 tmpgrad2 tmpgrad3 tmpgrad4 tmpgrad5 tmpgrad6 tmpgrad7 tmpgrad8 > dispgrad.dat
