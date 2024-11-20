# Time Evolving Block Decimation (TEBD)

We will now introduce our first MPS algorithm: Time Evolving Block Decimation. This algorithm allows us to simulate the time evolution of a quantum many-body system, where the state is represented by an MPS. The algorithm can also be used to simulate the imaginary (or complex) time evolution of a quantum system, which can be used to find ground states or to compute finite temperature correlation functions. We will focus on the real-time (unitary) time evolution in this course. We will consider the example of a global quench, to test our code against exact results. 

## Trotter decomposition

The TEBD algorithm is based on the Trotter decomposition of the time evolution operator. This allows us to approximate the evolution by a discrete sequence of local gates. The TEBD algorithm essentially tells us how to apply these gates to update our MPS.

The real-time evolution operator for a Hamiltonian $H$ is the unitary operator $U(t) = e^{-iHt}$. The Trotter (or Trotter-Suzuki) decomposition proceeds by two steps. First, we discretize time into small time steps $\Delta t$. The time evolution operator can then be written as $U(n\Delta t) = U(\Delta t)^n$. This first step is exact. The second step is to approximate the evolution operator $U(\Delta t)$ for a small time-step by a sequence of local unitary gates acting on pairs of sites (this can be made more general). Our Hamiltonian can be written as a sum of local terms $H = \sum_i h_i$, and we can approximate the evolution operator as

$$
U(\Delta t) = \prod_{n \; \text{even}} e^{-i h_n \Delta t} \prod_{n \; \text{odd}} e^{-i h_n \Delta t}  + O(\Delta t^2).
$$

This is the first-order Trotter decomposition, and is shown diagramatically in {numref}`fig:trotter`. It is also possible to use higher-order Trotter decompositions (and this is typically a good idea), but we will stick to first-order for simplicity. Therefore, we can approximate each time step by a sequence of local gates acting on pairs of sites of the form $e^{-i h_n \Delta t}$. We have arranged the product such that we apply all odd terms first, and then all even terms. This is because all the odd (even) terms commute with each other, but do not commute with the even (odd) terms. 

```{figure} images/trotter.png
---
name: fig:trotter
width: 80%
align: center
---

???
```

## Applying gates to MPS

The main routine of the TEBD algorithm is then to apply a two-site unitary to the MPS.

## Global quench test

```{figure} images/tebd_test.png
---
name: fig:extract_schmidt
width: 60%
align: center
---

???
```