#echo 'M=0.1'
QIAOYIN2=0
for QIAOYIN in {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.5,2.0,4.0,6.0,8.0,10.0}
do
QIAOYIN2=$((QIAOYIN2+1))
cp init.f90 ../init_$QIAOYIN2.f90
sed -i -e 's/M=2.0/M='$QIAOYIN'/g' init.f90
mm
mv main main__$QIAOYIN2
rm init.f90
mv ../init_$QIAOYIN2.f90 init.f90
echo 'M='$QIAOYIN 'No.'$QIAOYIN2
done
