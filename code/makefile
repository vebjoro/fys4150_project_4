
prob_4:
	g++ prob_4.cpp src/state.cpp -std=c++11 -I include -larmadillo  -o main
	./main

prob_56:
	g++ prob_56.cpp src/state.cpp -std=c++11 -I include -larmadillo  -o main
	./main

prob_8:
	g++ -c prob_8.cpp src/state.cpp -I include -fopenmp -std=c++11
	g++ prob_8.o state.o -larmadillo  -o main
	./main

prob_9:
	python3 prob_9.py

mac8:
	g++ prob_8.cpp src/state.cpp -Xpreprocessor -fopenmp -std=c++11 -I include -I /usr/local/Cellar/libomp/15.0.4/include -L /usr/local/Cellar/libomp/15.0.4/lib -larmadillo -lomp -o main.exe
	./main.exe

make all:
	make prob_4
	make prob_56
	make prob_8
	make prob_9

make all_mac:
	make prob_4
	make prob_56
	make mac8
	make prob_9
