[![DOI](https://zenodo.org/badge/666365140.svg)](https://zenodo.org/badge/latestdoi/666365140)

Disclaimer
----------

This repository was created for the sole purpose of archiving research data in a way that makes all numerical procedures and evaluations reproducible.
This means that the provided source codes do not meet the usual standards of a proper software release and that some advanced technical know-how is required to use them at all. 
I will not provide technical support beyond this README file, mostly because parts of the C++ code base (e.g. the entirety of `SRC-Reynolds`) were not written by me.



Info
----

Due to file size constraints, NO RAW DATA IS INCLUDED. 
Instead, the repository contains all necessary input files to repeat the simulations.
Once the corresponding scientific paper has been published, I will update the plot scripts and add the DOI here. 



Prerequisites
-------------

The contained Python scripts use well-known packages included in all Scientific Python Distributions, e.g. Anaconda ( https://www.anaconda.com/download ).
Apart from that, they make use of a few functions contained in my publically available contact mechanics utilities repository named `cmutils`:
https://github.com/sintharic/cmutils (make sure you also read the `README.md` there!)

The code in `SRC-GFMD` requires a working installation of FFTW3 and is compiled on the command line by typing `make`, generating the executable `contMech.exe`.
Using a package manager is recommended to install FFTW3, e.g. `sudo apt-get install libfftw3-dev` on Debian-based Linux distros. 
If you install FFTW in a non-standard location, you might need to change the corresponding path variables in the `Makefile`.

The code in `SRC-Reynolds` requires (an outdated version of) the `Hypre` package and is compiled on the command line by typing `./compileCurrent.sh`, generating the executable `currentCalc.x`.
Note that this is not a cross-platform solution and will almost certainly require some tinkering to work on your machine.



Usage
-----

You can have a look at all the plotted data inside the `analysis/Fig...` directories by opening the `.pdf` files or by plotting the `.txt` files yourself.
To reproduce the underlying raw data, you will have to rerun all simulations and post-process them yourself, the procedure of which is detailed in the following.

After installing `SRC-GFMD`, run the compiled `contMech.exe` in all sub-directories of `friction` that include a file named `params.in`.
Note that this represents a total of more than 250 simulations, which may be run in parallel (each one uses only a single CPU thread), but may still take a few days to finish.
Once all those simulations are done, run `prepare_flowX.py` and then `prepare_flowY.py` in the `analysis` directory.
Now you can run the Reynolds solver in all sub-directories of `flow` that contain a file named `gap_...`.
(To run the Reynolds solver in a folder, you have to copy the executable `currentCalc.x` from `SRC-Reynolds` to that folder and run `./currentCalc.x -n 2048 -solver 17 -iterMax 100000 > output.log &`.)
Once all of those simulations have finished, run all Python scripts in `analysis` whose file name starts with `res`. 
Then you can run all .py files in the `analysis/Fig...` folders to obtain the exact same figures as in the .pdf files.



Caveats
-------

Unfortunately, the behavior of the random number generator (RNG) in C++ is not standardized but implementation-specific. 
Therefore, even using the same RNG seed, the rough surface your machine generates might be different from the one that was generated for the publication.
Therefore, the original rough surface file `friction/rough/konfig0E.zip` is included.
