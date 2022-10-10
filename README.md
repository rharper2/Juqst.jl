# Juqst: JUst another Quantum Software Toolbox [![Build Status](https://travis-ci.com/rharper2/Juqst.jl.svg?branch=master)](https://travis-ci.com/rharper2/Juqst.jl)
# To install

It should be a matter of simply using the package manager.

This has not yet been put on the main public distribution, so you will need to point to the registry:

```julia
pkg> add https://github.com/rharper2/Juqst.jl
```

The exampleNotebooks directory contains, suprisingly enough, some sample notebooks that show how to work this. They assume IJulia. The one that shows how to use the Stabiliser mechanisms and plotting functions is called "A stabiliser run through". It is probably worth running that one early.


While you are able to view the notebooks on github, the latex rendering of printed arrays doesn't appear to work with the github viewer, so they might be a bit difficult to read on github. If you download them and look at them through Jupyter they should render fine. As mentioned above the workbook "A stabiliser run through" goes though the stabiliser mechanisms and plotting functions and is probably worth looking at if you are interested in the Stabiliser part of this package.

## Package Aims

The package consolidates a number of things I use when researching better methods of characterizing noise in quantum systems (often called QCVV for Quantum Control Verification and Validation). Just now it consists of three main branches:

- **Quantum Noise** - which contains all the algorithms and mechanisms necessary to extract $2^n$ eigenvalues/probabilities from a Quantum System in a single randomized benchmarking style experiment. It basically contains all the analysis code used in **Efficient Learning of quantum Noise** [arXiv:1907.13022](https://arxiv.org/abs/1907.13022). The workbooks, and graph drawing software all appear in **docs/example/quantumNoise** and show how to reproduce all the analysis and charts in that paper (and more), using the software in this package.
- **Channel conversions and random channel creation**. The first leans heavily on a package [open-system.jl]( https://github.com/BBN-Q/QuantumInfo.jl/blob/master/src/open-systems.jl) by Blake Johnson and Marcus da Silva - although rather than just import it, I have incorporated the software so I can tweak some minor aspects of it to fit in with the rest of my sofware. There are a number of minor convenience functions I have added, including the rendering I have added for arrays, the ability to generate 'Super-Vectors' and various labels for charts etc. The second part (in rchannels.jl) deals with ways of creating random noise channels, and utilised the insights contained in a paper by Rudnicki, Puchała, Zyczkowski [Gauge invariant information concerning quantum channels](https://arxiv.org/abs/1707.06926). This allows the creation of CPTP maps that have a nice spread of various characteristics such as fidelity and unitarity.
- **CHP and Clifford creation**. This implements the Stabiliser formalism of Aaronson and Gottesman and the and the ability to efficiently create/select arbitrary Cliffords by Koenig and Smolin. The Stabiliser formalism is prevalent in Quantum Computing and is going to become more relevant in QCVV as we move towards implementing error detection and correction codes. We can use this work to create the circuits we need to run. For instance by setting up a Clifford in the stabiliser state, we can automatically produce the Qiskit commands to create the Clifford.

I have still to incorporate the recent Bravyi and Maslov [Hadamard-free circuits expose the structure of the Clifford group](https://arxiv.org/abs/2003.09412) paper into this code. Futher minor todos are to check we have valid stabiliser state before we try and decompose (invalid stabiliser states can cause infinite loops) and to automate the finding of the destabilisers.

## Documentation 

This readme mainly details the CHP part of the programs included, there is documentation on the other functionality which can be found here: https://rharper2.github.io/Juqst.jl/docs/build/index.html. As perviously mentioned the docs/examples directory contains Jupyter notebooks that work through a lot of the functionality.

Below is a introduction to the CHP part of the package.

```julia
    Tableau(n::Integer)
```



 sets up and returns the a CHP Tableau for n qubits.

 This is based on the formalism given by: *Improved Simulation of Stabilizer Circuits*,
 Scott Aaronson and Daniel Gottesman, arXiv:quant-ph/0406196v5

 The initial tableau represents a |00...0⟩ ket in the stabiliser state
 This stabilises with "Z" and anti-stabilises with "X"

 For the purposes of this port, the tableau is exactly replicated as per the paper
 i.e. the "state" (Tableau.state) is an Int32 array (used as a bit array)
 containing the following information.

```
 x11   .....  x1n  | z11   ...       z1n | r1
  .    \\       .  |  .    \\          . |  .
  .     \\      .  |  .     \\         . |  .     Destabilisers
  .      \\     .  |  .      \\        . |  .
 xn1      \\   xnn | zn1      \\      znn| rn
 ______________________________________________
 x(n+1)1. x(n+1)n  | z(n+1) ... z(n+1)n  | r(n+1)
  .    \\       .  |  .      \\        . |  .
  .     \\      .  |  .       \\       . |  .     Stabilisers
  .      \\     .  |  .        \\      . |  .
 x(2n)1   \\x(2n)n | z(2n)1     \\ z(2n)n| r(2n)
```
Set Tableau.showRaw to true to see the underlying state as a matrix.

# Sample use

## Stabiliser Circuits

    state = setup(number_ofQubits)

prepares the stabiliser state for the correct number of qubits in the |000..000⟩ basis state

The state is represented internally as a matrix of the form:

<img src="readMeFigures/Matrix.png"></img>
Aaronson/Gottesman arXiv:quant-ph/0406196

Currently I am just using Int32 Arrays, although binary arrays would save space (if it ever becomes necessary).
Rows 1 to n of the tableau represent the destabiliser generators, rows n+1 to 2n represent the stabiliser generators. Each row is read
as follows: if the x<sub>ij</sub> and z<sub>ij</sub> are 1, the de/stabiliser is a Y, if they are both 0, its I otherwise its an X or Z depending on which one is set.

    output(state)

Prints the state in a human readable form. The states above the line are the 'destabiliser' state, below the line are the 'stabiliser' states. 

So in a 3 qubit system the initial state of |000> is coded as 

```
XII
IXI
IIX
---
ZII
IZI
IIZ
```

The following commands are defined

```JULIA
hadamard(t::Tableau,qubit)  # apply a hadamard to the relevant qubit
phase(t::Tableau,qubit)     # apply a phase gate to the relevant qubit
cnot(t::Tableau,control,target) # apply a controlled not from control qubit to target qubit
```

Output of the resultant state can be enabled by adding an extra true parameter

    hadamard(t::Tableau,qubit,true) # hadamard as before, but show output

**NOTE that these commands alter the state passed into them. I have broken Julia convention which requires functions 
with side effects to be written thus - hadamard!(state,qubit), however I think it is clear enough, it is after all the intended use of the function.**

## Arbitrary cliffords

(Koenig/Smolin arXiv:quant-ph/1406.2170)

The idea behind this paper is that we can implement a one-to-one mapping between the cliffords and an integer (plus a random phase string).

The mapping is as follows:

<img src="readMeFigures/Clifford Mapping.png">Koenig/Smolin arXiv:quant-ph/1406.2170</img>

We can generate the alpha,beta,gamma and delta via

   symplectic(i,n) # i = integer represting the clifford, n is the number of qubits

Which returns the nxn arrays (alpha->delta) coded as follows:

<img src="readMeFigures/coding.png">Koenig/Smolin arXiv:quant-ph/1406.2170</img>

More usefully these can be placed into a stabiliser tableau (that is the equivlent of passing the state |0000⟩ through a gate that implements the unitary in question as follows:

e.g.

    t = cliffordToTableau(4,23,1)

Where the qubits are 4 the Clifford chosen is 23, and we have chosen the first of $4^n$ phase patterns (here n = 4).

# Decomposing a tableau (such as clifford)

This will be made more general, but just now it decomposes an arbitrary tableau

    decompose(tableau)
    
    There is an optional parameter rationalise that defaults to true. Rationalise simply eliminates 4 phases in a row, two hadamards in a row or self cancelling cnots.

This prints out the elementary gates that would reconstruct the relevant clifford unitary.

The commands are stored as string in the vector commands


# Draw the circuit

This is a bit more involved, just now as it uses python packages. There should be an example in the notebooks. More details to be added.







