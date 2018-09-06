balance_Kaon_id : balance_Kaon_id.o
	g++ -O3 -m64 balance_Kaon_id.o -L/usr/local/root_v5.27/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lfreetype -pthread -Wl,-rpath,/usr/local/root_v5.27/lib -lm -ldl  -o balance_Kaon_id
balance_Kaon_id.o : balance_Kaon_id.C
	g++  -O3 -pipe -Wall -W -Woverloaded-virtual -D_REENTRANT -pthread -m64 -I/usr/local/root_v5.27/include -c balance_Kaon_id.C
