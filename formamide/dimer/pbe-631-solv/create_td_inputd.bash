#!/bin/bash

for i in {9..10}; do
    cat template.com geom.dat > td${i}.com
    sed -i "s/td1.chk/td${i}.chk/" td${i}.com
    sed -i "s/root=1/root=${i}/" td${i}.com
done
   
