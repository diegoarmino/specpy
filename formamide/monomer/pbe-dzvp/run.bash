#!/bin/bash

formchk freq.chk
for i in {1..6}; do
    g09 td${i}.com
    formchk td${i}.chk
    cp td${i}.fchk td.fchk
    ./calc_deltas.bin 6 0
    cp fort.77 disp${i}.dat
    grep -A12 gradient disp${i}.dat > tmpgrad${i}
done

paste tmpgrad* > dispgrad.dat
rm tmpgrad*
grep 'Excited State' td1.log > omega0.dat
grep -A12 'transition electric dipole moments' td1.log > trans-dip.dat
