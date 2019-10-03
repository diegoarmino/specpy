rm -f calc_deltas.bin
gfortran read_gaus_chk.f90 -L /usr/lib/x86_64-linux-gnu -lblas -o calc_deltas.bin
cp calc_deltas.bin gaus/

