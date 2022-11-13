comp:
	g++ main.cpp src/state.cpp -std=c++11 -I include/ -larmadillo -fopenmp -o main.exe

run:
	./main.exe

plt:
	python3 plot/plot_20x20.py

comprun:
	g++ main.cpp src/state.cpp -fopenmp -std=c++11 -I include/ -larmadillo -o main.exe
	./main.exe

comprunplt:
	g++ main.cpp src/state.cpp -std=c++11 -I include/ -larmadillo -fopenmp -o main.exe
	./main.exe
	python3 plot/plot_20x20.py

mac:
	g++ main.cpp src/state.cpp -Xpreprocessor -fopenmp -std=c++11 -I include/ -larmadillo -o -lomp main.exe
	./main.exe
