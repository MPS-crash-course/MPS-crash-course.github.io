# The Heisenberg Antiferromagnet

To guide our discussion of Matrix Product States (MPS) we will keep in mind a particular task: namely computing the dynamical structure factor for the Heisenberg Antiferromagnet. This task will combine all aspects of MPS that we will cover in this course. Let us start by introducing the model, provide some experimental motivation for this task, and then provide the details of the calculation we want to perform.

## The Heisenberg model

We will consider the one-dimensional Heisenberg model. This is a model of spin-1/2 particles (two-level quantum systems) on a one-dimensional chain. The spins are coupled to their nearest neighbours. The Hamiltonian for this model is given by

$$
H = J \sum_{i=1}^{N-1} \mathbf{S}_i \cdot \mathbf{S}_{i+1} = J \sum_{i=1}^{N-1} (S_i^x S_{i+1}^x + S_i^y S_{i+1}^y + S_i^z S_{i+1}^z)
$$ (eq:heisenberg)

where $S^\alpha_i = \frac{1}{2} \sigma^\alpha_i$ are the spin operators at site $i$, and $\sigma^\alpha_i$ are the Pauli matrices. $J$ is the coupling constant, and we will consider the case where $J=1$. In this case, the model is said to be antiferromagnetic since the spins prefer to be anti-aligned with their neighbours. The Heisenberg model is a paradigmatic model for quantum magnetism, which can be solved exactly in one dimension using Bethe ansatz. However, in this course we will use MPS to study the model. Crucially, the solvability can easily be broken by local perturbations, whereas the MPS techniques we will learn can be applied more generally.

## Experimental motivation

The Heisenberg model is an extremely simplified model for interacting spins, both because of the simplicity of the isotropic and uniform coupling, but also because it is a one-dimensional model. Remarkably, experimental measurements have shown that certain real materials are very accurately described by such a simple model. In the landmark work by B. Lake *et al* {cite}`Lake2013` the team measured the compound $\text{KCuF}_3$ using inelastic neutron scattering. By scattering high-energy neutrons off the sample, their deflection reveals information about the low energy excitations, as shown in Fig.~. The data is compared to the predictions of the antiferromagnetic Heisenberg model, and the agreement is remarkable. 

```{figure} ../../images/AFH_experiment.png
---
name: lake2013
width: 75%
align: center
---
Data from B. Lake *et al* {cite}`Lake2013`. The left panel shows the experimental data for the dynamical structure factor of $\text{KCuF}_3$. The right panel shows the theoretical prediction for the Heisenberg model obtained using Bethe Ansatz
```

This data reveals that the model is gapless, and has a very distinct dispersion relation. There is a single magnon (spin wave excitation) that has a dispersion that touches at $k=0, \pi$. Above this there is a continuum of mult-magnon excitations. The end goal of this course is to provide the theory predictions for this experiments using MPS methods.

## The dynamical structure factor



---

## References

```{bibliography}
:filter: docname in docnames
```