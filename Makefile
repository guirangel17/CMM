
movGen: communityMovGen-rel0.3.o socialNet.o
	g++ -o movGen communityMovGen-rel0.3.o socialNet.o

communityMovGen-rel0.3.o: communityMovGen-rel0.3.cc
	g++ -c communityMovGen-rel0.3.cc

socialNet.o: socialNet.cc
	g++ -c socialNet.cc
	
