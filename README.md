# Quantum-Decoherence

This repository contains several scripts that I have used to study the phenomenon of quantum decoherence; my results are a personal revisitation of previous results contained in important literature in the field.  
I presented them in the past, in [the slide](https://github.com/Edoardo1997/Quantum-Decoherence/blob/main/Presentation.pdf) (in Italian) I explain in detail the extent of my work.
# Brief explanation of my work

I have studied numerically the phenomenon of decoherence, one of the mechanisms that allow for the emergence of classical reality from the quantum world.  
The system under analysis is a harmonic oscillator in touch with a thermal bath. I want to show how the interaction with some external environment can force the selection of a preferred basis with classical properties. I solved the equation that rules this system (Caldeira-Legget master equation) numerically. 

The equation is a PDE in three variables and requires some attention to be solved maintaining numerical stability. I found out that the alternating direction Crank-Nicolson algorithm works perfectly and developed it. 
The algorithm computes the evolution of squeezed states showing how all states with squeezing die very fastly and a natural base of coherent states (which are famously states with classical properties) is selected.

<img src="/images/decoherence_filter.png" alt="drawing" width="400"/>

Naturally, the algorithm that I developed is more general and could be used to evolve any states, not just squeezed ones, but this was not of interest during my presentation.

It's interesting to know what happens after a natural basis is selected by the interaction with the environment. This doesn't require numerical calculation and simple analytical consideration do the trick; it turns out that when looking in the preferred base when you have quantum superpositions of states of the basis, they naturally lose their coherence terms (hence the term decoherence) and become simply statistical mixtures of classical states.

A visual representation of this phenomenon can be seen by plotting the density matrix of the quantum superposition of coherent states (the basis found before).

<p float="left">
  <img src="/images/cat_init.png" width="300" />
  <img src="/images/cat_final.png" width="300" /> 
</p>

The coherence terms vanish exponentially over time.  
Similarly, the same phenomenon can be noticed in the Wigner representation (where it's more clear what is classical and what is not as negative values are an indicator of quantumness)

<p float="left">
  <img src="/images/Wigner_init.png" width="300" />
  <img src="/images/Wigner_final.png" width="300" /> 
</p>
  
# How is this repository structured?

There are two main files that are relevant, everything else produces other material which was relevant for my presentation, but that is not related to the Crank-Nicolson algorithm.  
1. Decoherence.py: This file is responsible for the main computations, while real numerical simulations can take several hours, I've left some parameters which allows for a fast run (a few seconds). The data produced are stored in a file of byte with the Pickle library of Python, they can be plotted by View.py.
2. View.py: This file will plot the evaluated data. What it shows is a 3D graph which displays how certain states (the non-coherent ones) vanish faster than others (the coherent one) over time when they evolve according to Caldeira-Legget equations, thus showing how a base is naturally selected (decoherence).

There are other files responsible for the creation of pictures that I have used in my presentation, they only produce pictures that derive from analitycal calculations and not from the numerical solution of Caldeira-Legget master equation in Decoherence.py:
1. Cat_evol.py: Show the effect of quantum decoherence on the coherence terms of a quantum superposition that derives from the analytical solutions of the master equation in this case. 
2. Wigner.py: Show the same effect, but instead of using density matrix to represent the states, I used Wigner representation.

# How to run the code?
1. Clone this repository and move inside it
```bash
$ git clone https://github.com/Edoardo1997/Quantum-Decoherence.git
$ cd Quantum-Decoherence
```
2. (Suggested) Create a virtual environment
```bash
$ python3 -m venv quantumvenv
```
and activate it
```bash
$ source quantumvenv/bin/activate
```
3. Install requirements
```bash
$ pip install -r requirements.txt
```
4. To create new data: from terminal navigate to this directory and then run 
```bash
$ python3 source/Decoherence.py
```
5. To visualize the produced data: from this directory run
```bash
$ python3 source/View.py
```
I have also put some data from a more serious calculation (they took a couple of hours to be produced), if you want to visualize them, just change the path in View.py : /pickle_out.txt -> /save/pickle_out.txt.

6. To create the other pictures that I have used 
```bash
$ python3 images-source/Cat_evol.py
```
```bash
$ python3 images-source/Wigner.py
```
7. Exit from the virtual environment
```bash
$ deactivate
```
## References    
- Authors: Wojciech H. Zurek  
  Title: “Decoherence and the Transition from Quantum to Classical - Revisited”  
  In: Progress in Mathematical Physics, vol 48 (2006).  
  DOI: [10.1007/978-3-7643-7808-0_1](https://doi.org/10.1007/978-3-7643-7808-0_1)  

- Authors: Heinz-Peter Breuer and Francesco Petruccione  
  Title: “The Theory of Open Quantum Systems”  
  In: Oxford University Press (2002).  
  DOI: [0.1093/acprof:oso/9780199213900.001.0001](http://dx.doi.org/10.1093/acprof:oso/9780199213900.001.0001)
  
- Authors: Frank Grossmann and Werner Koch  
  Title: “A finite-difference implementation of the Caldeira–Leggett master equation”    
  In: J. Chem. Phys. 130, 034105 (2009).  
  DOI: [10.1063/1.3059006](https://doi.org/10.1063/1.3059006)    
       
