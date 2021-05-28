# Quantum-Decoherence

This repository contains several script that I have used to study the phenomenon of quantum decoherence, my results are a personal revisitation of previous results contained in important literature in the field. I exposed in the past my results in a presentation [the presentations slide](https://github.com/Edoardo1997/Quantum-Decoherence/blob/main/Presentation.pdf) (in Italian) explains in detail the extent of my work
# What have I done?

I have developed an alternating direction Crank-Nicolson algoirthm to solve numerically Caldeira-Legget equations to study the evolution of an harmonic oscillator in touch with an external bath. I showed numerically that when the system is underdamped the environment act naturally selecting a preferred base of the Hilbert space (more specifically the set of coherent states which are overcomplete).
# What is in the various file?

There are two main files which are relevant, everything else produce other material which was relevant for my presentation but that is not related to the Crank-Nicolson algorithm.
1. Decoherence.py : This file is responsible for the main computations, while real numerical simulations can take several hours, I've left some parameters which allows for a fast run (a few seconds). The data produced are stored in a files of byte with the Pickle library of Python and they can be plotted by view.py.
2. view.py : This file will plot the evaluated data, what is shown is a 3D graph which shows how certain states (the non-coherent ones) goes to zero faster than others (the coherent one) over time when they evolve according to Caldeira-Legget equations, thus showing how a base is naturally selected (decoherence).

# How to run my script?
1. To create new data: from terminal navigate to this directory and then run 'python3 Decoherence.py'.
2. To visualize the produced data: from this directory run 'python3 view.py'.
I have put also some data from a more serious calculation (they took a a couple of hours to be produced), if you want to visualize them just change the path in 'view.py' : '/pickle_out.txt' -> '/save/pickle_out.txt'. 
