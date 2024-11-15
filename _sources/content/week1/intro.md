# Introduction

Before getting into what we will be doing in the course, let me start with a bit of context and demystify some of the jargon.

## Tensor Networks (TNs)

Before getting to Matrix Product States, let's start at the highest level: Tensor Networks (TNs). Tensor networks are an extremely general way of representing different objects and processes, and are used in many areas of physics, chemistry, and mathematics. In the most abstract sense, they are represented by a graph. Each vertex of the graph is a tensor, and each edge is a tensor contraction. This will hopefully become clearer when we deal with the details of Matrix Product States. 

In physics, they are most commonly used as a way to represent (or approximate) quantum many-body states. However, they can also be used to represent classical probability distributions, quantum channels, and many other things. Quantum circuits are also naturally represented as a type of TNs. 

Matrix Product States are a specific type of one-dimensional tensor network. However, there are also Projected Entangled Pair States (PEPS), Multi-scale Entanglement Renormalization Ansatz (MERA), Tree Tensor Networks (TTN), and many others.

If you are interested in matrix product states and tensor networks more generally then I would recommend the following reviews {cite}`Schollwock2011,Orus2014,Cirac2021`.

## Matrix Product States (MPS)

Matrix Product States (MPS) are a particular type of TN that has an intrinsic one-dimensional structure and are the focus of this course. They are used to represent quantum many-body states in one dimension, but can also be used in higher dimensions.

The key idea behind MPS (and indeed the origin of the name) is that the probability amplitudes of a quantum state can be represented as a product of matrices. Importantly, this representation allows us to truncate the size of the matrices, which allows us to controllably *approximate* quantum states, making computations tractable. Exactly how we should truncate the matrices is intrinsically related to entanglement entropy, and we will discuss this in week 2. 

The modern way to view MPS is as a variational ansatz for quantum many-body states, where the size of the matrices is restricted. This allows us to study very large systems (>100 sites) that would be intractable with exact diagonalization, at the cost of restricting to a subspace of the full Hilbert space, and therefore potentially inducing errors if the desired states are not well approximated by the MPS ansatz. The power of MPS is that they are provably good at efficiently representing ground states of one-dimensional quantum systems.

Matrix product states were first introduced in the context of the Affleck-Kennedy-Lieb-Tasaki (AKLT) model {cite}`Affleck1988`, where the exact ground state can be represented as a product of finite size matrices {cite}`Klumper1992`. 


## Density Matrix Renormalization Group (DMRG)

Alongside the exact representation of the AKLT state, the second main ingredient that led to modern MPS methods is the Density Matrix Renormalization Group (DMRG) algorithm {cite}`White1992`. DMRG is a numerical method for finding the ground state of one-dimensional quantum systems. The original paper did not use the language of matrix product states, but instead built up a many-body state using a real-space renormalization group approach. The key insight came in how to truncate the Hilbert space in a controlled way. 

```{figure} images/dmrg.png
---
name: fig:dmrg
width: 40%
align: center
---

Schematic of the original DMRG algorithm, taken from {cite}`White1992`. The system is built up by adding two sites at a time, and the reduced density matrix of the half-chain is used to truncate the Hilbert space.
```

The main idea was to successively increase the system size by adding two sites to the centre of the system. By looking at the reduced density matrix for the half chain, only the states with the largest probability are kept, as shown in {numref}`fig:dmrg`. This truncation maximises the allowed entanglement between the two halves of the system. The process can then be repeated many times since only a finite number of states is ever kept, allowing you to approximate the ground state energy (and state) for large or infinite chains. This truncation in terms of the entanglement spectrum underlies all MPS methods. The DMRG in now used in a form expressed in terms of MPS, and is one of the most successful numerical methods for studying ground states of one-dimensional quantum systems. **We will discuss DMRG in week 3**.


## Time Evolving Block Decimation (TEBD)

MPS are most powerful for studying ground state of one-dimensional quantum systems due to their particular entanglement structure (a topic we won't go into detail about). However, it was also realized that they can be used to study the dynamics of quantum systems, and in many cases allow us to get useful results that are not accessible using exact diagonalization based methods.

While there are several ways of performing time evolution, the one we will consider, and the most historically important, is the Time Evolving Block Decimation (TEBD) algorithm {cite}`Vidal2003,Vidal2004`. This is a method for evolving a matrix product state in time. The key idea is to represent the time evolution operator as a series of local gates, which can be applied to the MPS. This allows us to simulate the time evolution of a quantum state. During the evolution we are able to perform the same truncation as in DMRG, which allows us to simulate the dynamics of large quantum systems. By evolving using imaginary time, TEBD can also be used to find ground states, or to compute thermal expectation values. However, DMRG is far more efficient for finding ground states in practice, and TEBD is mainly used for dynamics. **We will discuss TEBD in week 2**.


---

## References

```{bibliography}
:filter: docname in docnames
```


