formchk freq.chk
for i in {1..2}; do 
    formchk td${i}.chk
    cp td${i}.fchk td.fchk
    #cp td${i}-opt.fchk td-opt.fchk
    ./calc_deltas.bin 6 0
    cp fort.77 disp${i}.dat
done

for i in {1..2}; do 
    #grep -A12 LMXg disp${i}.dat > tmpopt${i}
    grep -A12 gradient disp${i}.dat > tmpgrad${i}
done

#paste tmpopt1 tmpopt2 tmpopt3 tmpopt4 tmpopt5 tmpopt6 tmpopt7 tmpopt8 > dispopt.dat
paste tmpgrad1 tmpgrad2 > dispgrad.dat
