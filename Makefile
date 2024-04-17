all:
	g++ -march=native -O3 -std=c++11 -o triangle_DP triangle_DP.cpp
	g++ -march=native -O3 -std=c++11 -o triangle_LDP triangle_LDP.cpp

clean:
	rm -f triangle_DP triangle_LDP
