# The Heisenberg Antiferromagnet

To guide our discussion of Matrix Product States (MPS) we will keep in mind a particular task: namely computing the dynamical structure factor for the Heisenberg Antiferromagnet. This task will combine all aspects of MPS that we will cover in this course. Let us start by introducing the model, provide some experimental motivation for this task, and then provide the details of the calculation we want to perform.

## The Heisenberg model

We will consider the one-dimensional Heisenberg model. This is a model of spin-1/2 particles (two-level quantum systems) on a one-dimensional chain. The spins are coupled to their nearest neighbours. The Hamiltonian for this model is given by

$$
H = J \sum_{i=1}^{N-1} \mathbf{S}_i \cdot \mathbf{S}_{i+1} = J \sum_{i=1}^{N-1} (S_i^x S_{i+1}^x + S_i^y S_{i+1}^y + S_i^z S_{i+1}^z)
$$ (eq:heisenberg)

where $S^\alpha_i = \frac{1}{2} \sigma^\alpha_i$ are the spin operators at site $i$, and $\sigma^\alpha_i$ are the Pauli matrices. $J$ is the coupling constant, and we will consider the case where $J=1$. In this case, the model is said to be antiferromagnetic since the spins prefer to be anti-aligned with their neighbours. The Heisenberg model is a paradigmatic model for quantum magnetism, which can be solved exactly in one dimension using Bethe ansatz. However, in this course we will use MPS to study the model. Crucially, the solvability can easily be broken by local perturbations, whereas the MPS techniques we will learn can be applied more generally.

## Experimental motivation

The Heisenberg model is an extremely simplified model for interacting spins, both because of the simplicity of the isotropic and uniform coupling, but also because it is a one-dimensional model. Remarkably, experimental measurements have shown that certain real materials are very accurately described by such a simple model. In the landmark work by B. Lake *et al* {cite}`Lake2013` the team measured the compound $\text{KCuF}_3$ using inelastic neutron scattering. By scattering high-energy neutrons off the sample, their deflection reveals information about the low energy excitations, as shown in {numref}`fig:lake2013`. The data is compared to the predictions of the antiferromagnetic Heisenberg model, and the agreement is remarkable. 

```{figure} images/AFH_experiment.png
---
name: fig:lake2013
width: 80%
align: center
---
Data from B. Lake *et al* {cite}`Lake2013`. The left panel shows the experimental data for the dynamical structure factor of $\text{KCuF}_3$. The right panel shows the theoretical prediction for the Heisenberg model obtained using Bethe Ansatz.
```

This data reveals that the model is gapless, and has a very distinct dispersion relation. There is a single magnon (spin wave excitation) that has a dispersion that touches at $k=0, \pi$. Above this there is a continuum of multi-magnon excitations. The end goal of this course is to provide the theory predictions for this experiments using MPS methods.

## The dynamical structure factor

The specific quantity that we want to compute is called the *dynamical spin structure factor (DSSF)*. This is the quantity that can be measured in inelastic neutron scattering experiments. The components of the DSSF are given by

$$
S^{\alpha,\beta}(q, \omega) = \frac{1}{2\pi N}\sum_{j,k} e^{iq(j-k)} \int_{-\infty}^{\infty} e^{i\omega t} \langle S^\alpha_j(t) S^\beta_k(0) \rangle dt.
$$ (eq:dsf)

The actual quantity that is measured in experiments is the cross-section for the outgoing neutrons, which is related to the dynamical structure factor but also includes simple momentum polarization and magnetic form factors. These details go beyond the scope of this course. In fact, due to symmetries of the model, we will only be concerned with the $S^{zz}(q, \omega)$ component.

The DSSF is the Fourier transform of a two-point unequal time spin correlation function,

$$
C^{zz}(r, t) = \langle S^z_{j+r}(t) S^z_{j}(0) \rangle = \langle \psi_0 | e^{iHt} S^z_{j+r} e^{-iHt} S^z_{j} | \psi_0 \rangle
$$ (eq:corr)

where $|\psi_0\rangle$ is the ground state of the Antiferromagnetic Heisenberg model, described by the Hamiltonian $H$. This quantity is in general prohibitively costly to compute with exact numerics due to the exponential growth of the Hilbert space for quantum states, and is an ideal candidate for MPS methods.
While we will go into the details of this calculation in week 4, our goal can then be split into two main parts: 
* computing the ground state $|\psi_0\rangle$ of the Heisenberg model. We will do this using the Density Matrix Renormalization Group (DMRG) method. 
* performing the unitary time evolution under $e^{-iHt}$, which will be achieved by using the Time-Evolving Block Decimation (TEBD) algorithm.


By the end of the course you should be able to produce data similar to that shown in {numref}`fig:mpsfinal` using your own code:

```{figure} images/dssf_100_16.jpg
---
name: fig:mpsfinal
width: 66%
align: center
---
Data for the dynamical spin structure factor of the Heisenberg Antiferromagnet obtained using MPS methods you will learn in this course.
```


````{admonition} Code: Basic Exact Diagonalization

Throughout this course, it will be good to check our MPS results against exact diagonalization. Since I do not expect you to have done exact diagonalization before, I will provide the code. For now, the code will just allow us to construct the Hamiltonian matrix for the Heisenberg model and to compute the ground state energy and wavefunction. We will extend this throughout the course as necessary.

```python
## file: src/ed.py

import numpy as np
import scipy.sparse as sp


def HeisenbergHamiltonian(L):
    """
    Construct the Heisenberg Hamiltonian for a 1D chain of length L.
    """
    # Define the spin operators
    s_x = 1/2*np.array([[0, 1], [1, 0]])
    s_y = 1/2*np.array([[0, -1j], [1j, 0]])
    s_z = 1/2*np.array([[1, 0], [0, -1]])

    # Construct the Heisenberg Hamiltonian
    H = np.zeros((2**L, 2**L))
    for i in range(L-1):
        H += np.kron(np.kron(np.kron(np.eye(2**i), s_x), s_x), np.eye(2**(L-i-2)))
        H += np.real(np.kron(np.kron(np.kron(np.eye(2**i), s_y), s_y), np.eye(2**(L-i-2))))
        H += np.kron(np.kron(np.kron(np.eye(2**i), s_z), s_z), np.eye(2**(L-i-2)))

    return H


def HeisenbergGroundState(L):
    """
    Compute the ground state of the Heisenberg Hamiltonian for a 1D chain of length L.
    """
    H = HeisenbergHamiltonian(L)
    E, V = sp.linalg.eigsh(H, k=1, which='SA')
    
    return E[0], V[:,0]

```

````


---

## References

```{bibliography}
:filter: docname in docnames
```