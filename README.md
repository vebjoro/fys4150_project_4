# Monte Carlo simulation of Ising model

This project aims to simulate the Ising Model numerically using Markov Chain Monte Carlo methods.

## How to run this code?
Make sure all dependencies are installed. To compile C++ files, run the following command:

#### MacOS:
```
make all_mac
```

#### Linux/Windows
```
make all
```
NOTE: These programs have only been tested on MacOS.

To run the pyton scrips, run the following commands with \/code` as working directory:
```
python3 plot plot/plot_20x20.py
```

```
python3 plot plot/plot_OMP.py
```


### Dependencies
- Armadillo (C++)
- OpenMP (C++)
- Matplotlib (Python)
- Numpy (Pyhton)
- Pyarma (Python)
- Scipy (Python)