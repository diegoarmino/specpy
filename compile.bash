rm -f calc_deltas.bin
gfortran calc_deltas.f90 -L/usr/lib/x86_64-linux-gnu/lapack -llapack -lblas -o calc_deltas.bin
cp calc_deltas.bin formamide/cam-b3lyp-tz-solv/

