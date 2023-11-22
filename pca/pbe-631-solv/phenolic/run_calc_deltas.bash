for i in {1..9}
do
	cp td${i}.fchk td.fchk
	./calc_deltas.bin 6 0
	grep -A12 gradient fort.77 > disp${i}.dat
done

