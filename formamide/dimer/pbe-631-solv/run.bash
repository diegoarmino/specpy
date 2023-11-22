#!/bin/bash

nstates=10
nat=12
ndf=36

for i in {2..10}; do
    g09 td${i}.com
done

formchk freq.chk
for i in {1..10}; do
    formchk td${i}.chk
    cp td${i}.fchk td.fchk
    ./calc_deltas.bin ${nat} 0
    cp fort.77 disp${i}.dat
    grep -A${ndf} gradient disp${i}.dat > tmpgrad${i}
done

paste tmpgrad* > tmpgradall
awk '{printf "%10s%10s%10s%10s%10s%10s\n",$1,$3,$5,$7,$9,$11,$13,$15,$17,$19,$20}' tmpgradall > dispgrad.dat
#rm tmpgrad*
#grep 'Excited State' td1.log > omega0.dat
#grep -A12 'transition electric dipole moments' td1.log > trans-dip.dat
