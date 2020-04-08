# Simplified Plotting Routines for Quantum Circuits
Rick Muller

This program takes some of the circuit tricks that Brian Granger, Aaron Meurer and Ondrej Certik first developed in Sympy and [I subsequently updated](http://nbviewer.jupyter.org/gist/rpmuller/5843312).

I pulled out all of the code into a set of standalone python functions because I wanted to:
* have fewer dependencies for people to install;
* make a simpler code base to experiment with;
* separate the quantum circuit from the plotter from the simulator;
* experiment with scheduling quantum operations rather than plotting one gate in each time step.

There are lots of other tools you can do for circuit plotting.
* [Qasm2circ](https://www.media.mit.edu/quanta/qasm2circ/)
* [Q-circuit](http://physics.unm.edu/CQuIC/Qcircuit/)
* The aforementioned [Sympy circuit plotter](http://nbviewer.jupyter.org/gist/rpmuller/5843312).

## TODO List:
* Figure out an elegant way to put multiqubit targets into the gates.
* Use multiqubit rectangles for multitarget gates, if the wires are adjacent.

