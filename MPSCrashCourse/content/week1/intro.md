# Introduction

Before getting into what we will be doing in the course, let me start with a bit of context and demystify some of the jargon.

## Tensor Networks (TNs)

Before getting to Matrix Product States, let's start at the highest level: Tensor Networks (TNs). Tensor networks are an extremely general way of representing different objects and processes, and are used in many areas of physics, chemistry, and mathematics. In the most abstract sense, they are represented by a graph. Each vertex of the graph is a tensor, and each edge is a tensor contraction. This will hopefully become clearer when we deal with the details of Matrix Product States. 

In physics, they are most commonly used as a way to represent (or approximate) quantum many-body states. However, they can also be used to represent classical probability distributions, quantum channels, and many other things. Quantum circuits are also naturally represented as a type of TNs. 

Matrix Product States are a specific type of one-dimensional tensor network. However, there are also Projected Entangled Pair States (PEPS), Multi-scale Entanglement Renormalization Ansatz (MERA), Tree Tensor Networks (TTN), and many others.

...Add references for tensor network reviews!...

## Matrix Product States (MPS)

Matrix Product States (MPS) are a particular type of TN that has an intrinsic one-dimensional structure and are the focus of this course. They are used to represent quantum many-body states in one dimension, but can also be used in higher dimensions.

The key idea behind MPS (and indeed the origin of the name) is that the probability amplitudes of a quantum state can be represented as a product of matrices. Importantly, this representation allows us to truncate the size of the matrices, which allows us to controllably *approximate* quantum states, making computations tractable. Exactly how we should truncate the matrices is intrinsically related to entanglement entropy, and we will discuss this in week 2. 

The modern way to view MPS is as a variational ansatz for quantum many-body states, where the size of the matrices is restricted. This allows us to study very large systems (>100 sites) that would be intractable with exact diagonalization, at the cost of restricting to a subspace of the full Hilbert space, and therefore potentially inducing errors if the desired states are not well approximated by the MPS ansatz. The power of MPS is that they are provably good at efficiently representing ground states of one-dimensional quantum systems.

Matrix product states were first introduced in the context of the Affleck-Kennedy-Lieb-Tasaki (AKLT) model {cite}`Affleck1988`, where the exact ground state can be represented as a product of finite size matrices {cite}`Klumper1992`. 


## Density Matrix Renormalization Group (DMRG)




---

## References

```{bibliography}
:filter: docname in docnames
```


