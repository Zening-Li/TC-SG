all:
	g++ -march=native -O3 -std=c++11 -o DP DP.cpp
	g++ -march=native -O3 -std=c++11 -o DPSaveMem DPSaveMem.cpp
	g++ -march=native -O3 -std=c++11 -o LDP LDP.cpp

clean:
	rm -f DP DPSaveMem LDP
