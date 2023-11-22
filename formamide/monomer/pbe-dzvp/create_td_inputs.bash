for i in {1..8};
do 
   cat td.com.template geom.dat dzvp.gbs > td${i}.com
   sed -i "s/td1/td${i}/" td${i}.com
   sed -i "s/root\=1/root\=${i}/" td${i}.com
   sed -i 's/6-31G\*\*/Gen/' td${i}.com
   sed -i 's/cam-b3lyp/pbepbe/' td${i}.com
done
