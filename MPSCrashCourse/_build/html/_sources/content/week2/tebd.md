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




```{figure} images/tebd_test.png
---
name: fig:extract_schmidt
width: 60%
align: center
---

???
```



````{admonition} Exercise: Local Quench

A

```{figure} images/tebd_exercise.png
---
name: fig:tebd_exercise
width: 60%
align: center
---

???
```



````
