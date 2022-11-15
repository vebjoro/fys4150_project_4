
prob_4:
	g++ prob_4.cpp src/state.cpp -std=c++11 -I include -larmadillo  -o main
	./main

prob_56:
	g++ prob_56.cpp src/state.cpp -std=c++11 -I include -larmadillo  -o main
	./main

prob_8:
	g++ prob_8.cpp src/state.cpp -fopenmp -std=c++11 -I include -larmadillo  -lomp -o main
	./main

mac8:
	g++ prob_8.cpp src/state.cpp -Xpreprocessor fopenmp -std=c++11 -I include -I /usr/local/Cellar/libomp/15.0.4/include -L /usr/local/Cellar/libomp/15.0.4/lib -larmadillo -lomp -o main.exe
	./main.exe
