# Juqst Documentation

This is some software that brings together:
- the software and data that was used to analyse and craeate the plots for [Efficient Learning of Quantum Noise](https://www.nature.com/articles/s41567-020-0992-8) and [arXiv](https://arxiv.org/abs/1907.13022) (see examples/noise) 
- a simulation of stabiliser circuits [Aaronson/Gottesman arXiv:quant-ph/0406196](http://arxiv.org/pdf/quant-ph/0406196)
- code to select an arbitrary Clifford group element [Koenig/Smolin arXiv:quant-ph/1406.2170](http://arxiv.org/abs/1406.2170)
- decomposition an arbitrary clifford unitary into a quantum circuit consistiting of hadamard, phase and two-qubit cnot gates [Aaronson/Gottesman arXiv:quant-ph/0406196](http://arxiv.org/pdf/quant-ph/0406196)
- conversion of the Tableau to a unitary matrix
- a slightly modified version of Marcus da Silva's QuantumInfo.jl - updated with some minor modifications and additions.
- some additional functionality to allow the creation of random quantum channels (for testing purposes).


```@contents
Pages = ["probability-algs.md","chp.md", "open-systems.md","channels.md"]
Depth = 2
```
GitHub repo [Juqst.jl](https://github.com/rharper2/Juqst.jl)
