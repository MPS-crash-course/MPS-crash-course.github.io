# Time Evolving Block Decimation (TEBD)

We will now introduce our first MPS algorithm: Time Evolving Block Decimation. This algorithm allows us to simulate the time evolution of a quantum many-body system, where the state is represented by an MPS. The algorithm can also be used to simulate the imaginary (or complex) time evolution of a quantum system, which can be used to find ground states or to compute finite temperature correlation functions. We will focus on the real-time (unitary) time evolution in this course. We will consider the example of a global quench, to test our code against exact results. 

## Trotter decomposition

The TEBD algorithm is based on the Trotter decomposition of the time evolution operator. This allows us to approximate the evolution by a discrete sequence of local gates. The TEBD algorithm essentially tells us how to apply these gates to update our MPS.

The real-time evolution operator for a Hamiltonian $H$ is the unitary operator $U(t) = e^{-iHt}$. The Trotter (or Trotter-Suzuki) decomposition proceeds by two steps. First, we discretize time into small time steps $\Delta t$. The time evolution operator can then be written as $U(n\Delta t) = U(\Delta t)^n$. This first step is exact. The second step is to approximate the evolution operator $U(\Delta t)$ for a small time-step by a sequence of local unitary gates acting on pairs of sites (this can be made more general). Our Hamiltonian can be written as a sum of local terms $H = \sum_i h_i$, and we can approximate the evolution operator as

$$
U(\Delta t) = \prod_{n \; \text{even}} e^{-i h_n \Delta t} \prod_{n \; \text{odd}} e^{-i h_n \Delta t}  + O(\Delta t^2).
$$

This is the first-order Trotter decomposition, and is shown diagramatically in {numref}`fig:trotter`. It is also possible to use higher-order Trotter decompositions (and this is typically a good idea), but we will stick to first-order for simplicity. Therefore, we can approximate each time step by a sequence of local gates acting on pairs of sites of the form $e^{-i h_n \Delta t}$. We have arranged the product such that we apply all odd terms first, and then all even terms. This is because all the odd (even) terms commute with each other, but do not commute with the even (odd) terms. 

```{figure} images/trotter.jpeg
---
name: fig:trotter
width: 60%
align: center
---

???
```

## Applying gates to MPS

The main routine of the TEBD algorithm is then to apply a two-site unitary to the MPS as shown in {numref}`fig:theta`. First we must move the centre to one of the two site we are acting on. We will do sweeps from left to right, so we will move the centre to the left site. We then reshape the $4 \times 4$ unitary into a tensor with shape $(2,2,2,2)$. We can then contract this tensor with the two tensors it acts on (shown in left of {numref}`fig:theta`) to get a new tensor, which we call $\Theta$. Finally, we perform a truncated SVD on $\Theta$ to get the state back in MPS form with updated tensors. All that is left to complete one TEBD step is loop over all the odd unitary gates, then all the even gates.


```{figure} images/theta.jpeg
---
name: fig:theta
width: 100%
align: center
---

???
```

```{note}
Just a note on keeping track of index order in python. If we have a $4 \times 4$ unitary, I would label the indices (out, in). Meaning that I would contract the `in` index with a vector if I were doing matrix-vector multiplication. If I reshape this matrix into a rank-4 tensor then I would get (out_left, out_right, in_left, in_right), where out_left refers to the leg coming out of the tensor corresponding to the left spin. In {numref}`fig:theta` the in indices are at the top, and out at the bottom. Similarly for the other indices. WARNING: if you are using Julia, then the labelling would be different.

```


````{admonition} Code: TEBD Step

Let us create a new file for our TEBD code called `tebd.py` in the `src` folder. We will need to add the following functions. I leave the implementation of the `applyGate` function to you. 

```python
## file: src/tebd.py

import numpy as np
import scipy as sp

from .svd import svd_truncated

PauliX = np.array([[0,1],[1,0]])
PauliY = np.array([[0,-1j],[1j,0]])
PauliZ = np.array([[1,0],[0,-1]])
XX = np.kron(PauliX, PauliX)
YY = np.kron(PauliY, PauliY)
ZZ = np.kron(PauliZ, PauliZ)
heisenberg_term = 1/4*(XX + YY + ZZ)

def heisenbergGate(dt):
    """
    Compute the two-site Heisenberg gate for a given time step dt.
    """

    return sp.linalg.expm(-1j*dt*heisenberg_term)


def applyGate(psi, site, dt, chiMax, tol):
    """
    Apply a two-site gate for the Heisenberg model on site and site+1.
    """

    gate = heisenbergGate(dt)

    psi.move_centre_to(site)

    ## YOUR CODE HERE ##



def TEBD_step(psi, dt, chiMax, tol):
    """
    Perform a single time step of the time-evolution block decimation (TEBD) algorithm. First-order trotterization of the Heisenberg Hamiltonian.

    Parameters
    ----------
    psi : MPS
        Matrix Product State object.
    dt : float
        Time step.
    """

    L = psi.L

    # apply the Heisenberg gate to each pair of odd sites
    for i in range(0,L-1,2):
        applyGate(psi, i, dt, chiMax=chiMax, tol=tol)

    for i in range(1,L-1,2):
        applyGate(psi, i, dt, chiMax=chiMax, tol=tol)


```

````


## Global quench test

We now have all the ingredients to test our TEBD algorithm. Let us consider a global quench of the Heisenberg model. Starting from a product state with alternating spins up and down, i.e. $|0101010\cdots \rangle$, we will evolve this state under the Heisenberg Hamiltonian. We will measure the magnetization of the system at the central site as a function of time, as well as the entanglement entropy, as shown in {numref}`fig:tebd_test`. We will compare our results for different bond dimensions $\chi$ to the results from exact diagonalization.

```{figure} images/tebd_test.png
---
name: fig:tebd_test
width: 60%
align: center
---

???
```

In {numref}`fig:tebd_test` the simulation was performed for $L=10, dt=0.05, t_\text{max}=6$. The accuracy tolerance was set to `None` and the MPS code was run for different values of $\chi_\text{max} = 2,4,8,16$. We can see that the TEBD simulation matches the exact results well up to time set by $\chi_\text{max}$. The entanglement entropy shows that the point where the results deviate corresponds to the point where the entanglement saturates due to the SVD truncation. At this point, we are throwing too much information about the state away and the agreement gets worse. By running simulations with varying $\chi_\text{max}$, and by monitoring the entanglement entropy, we can determine up to which point we can trust our results.


````{admonition} ED Code: Global Quench

As before, I provide here the ED code for the global quench. This code will be used to compare the results of the TEBD algorithm. 

```python
## file: src/ed.py

## PREVIOUS CODE OMITTED ##

def entanglementEntropy(psi, site):
    """
    Compute the entanglement entropy of a quantum state psi across bond between site and site+1.
    """

    psi = psi.copy().reshape((2**(site+1), -1))
    _, S, _ = np.linalg.svd(psi, full_matrices=False)

    return -np.sum(S**2 * np.log(S**2))


def HeisenbergTimeEvolution(L, state, dt, tMax):
    """
    Compute the time evolution of the Heisenberg Hamiltonian for a 1D chain of length L. L is even!
    """
    H = HeisenbergHamiltonian(L)
    psi0 = 1.
    for i in state:
        psi0 = np.kron(psi0, np.array([int(i == 0), int(i == 1)]))
    
    psi = psi0
    nSteps = int(tMax / dt)

    U = expm(-1j*dt*H)
    Z = np.kron(np.kron(np.eye(2**(L//2)), np.array([[1,0],[0,-1]])),np.eye(2**(L//2-1)))

    magnetization = []
    entanglement = []
    for i in range(nSteps):
        psi = U @ psi
        magnetization.append( np.real( psi.conj().T @ Z @ psi ) )
        entanglement.append( entanglementEntropy(psi, L//2-1) )
    
    return magnetization, entanglement

```

````




````{admonition} Exercise: Local Quench

A good exercise for you would be to implement the TEBD algorithm for a local quench. Consider a 1D chain in an initial product state where all spins are up, except for three. I have chosen to place these at L/4, L/2, 3L/4. You can then perform the time evolution using TEBD and measure the magnetization on all sites, allowing you to plot the magnetization profile as a function of time, as shown in {numref}`fig:tebd_exercise`. 

```{figure} images/tebd_exercise.png
---
name: fig:tebd_exercise
width: 60%
align: center
---

???
```

To produce this figure I set $L=51, dt=0.1, t_\text{max}=25$. The accuracy tolerance was set to `None` and the MPS code was run for $\chi_\text{max} = 8$ (which gives exact results in this case). This showcases the ability for MPS methods to simulate time evolution for large systems sizes, without the need to store the full state vector.



````
