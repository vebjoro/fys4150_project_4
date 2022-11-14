linux:
	g++ main.cpp src/state.cpp  -fopenmp -std=c++11 -I include -larmadillo -lomp -o main.exe
	./main.exe

mac:
	g++ main.cpp src/state.cpp -Xpreprocessor -fopenmp -std=c++11 -I include -I /usr/local/Cellar/libomp/15.0.4/include -L /usr/local/Cellar/libomp/15.0.4/lib -larmadillo -lomp -o main.exe
	./main.exe