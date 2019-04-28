for k in 0.0 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 12.0 13.0 14.0 15.0 16.0 17.0 18.0 19.0 20.0; do
	 movGen -t 6000 -n 50 -s $k -S $k -X 1000.0 -Y 1000.0 -R 5 -C 5 -G 5 -w 0.1 -g 5 -a off -d off
         cp socMov.xml socMovCom-$k.xml 
done
