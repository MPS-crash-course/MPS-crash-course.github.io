# Bringing It All Together

At this point we have learnt all the aspects of MPS necessary to tackle the problem we set out at the start: computing the dynamical spin structure factor for the antiferromagnetic Heisenberg model. We have learnt about the MPS representation of quantum states, the DMRG algorithm for finding ground states, and the TEBD for simulating time evolution. It is now time to bring all these elements together.

## Dynamical Spin Structure Factor

Let us restate the problem that we are tackling with our MPS code. We want to compute the *dynamical spin structure factor* (DSSF)

$$
S^{\alpha\beta}(q, \omega) = \frac{1}{2\pi N}\sum_{j,k} e^{iq(j-k)} \int_{-\infty}^{\infty} e^{i\omega t} \langle \psi_0 | S^\alpha_j(t) S^\beta_k(0) | \psi_0 \rangle \;dt,
$$ (eq:dsf_2)

for the antiferromagnetic Heisenberg (AFH) model described by the Hamiltonian

$$
H = J \sum_{i=1}^{N-1} \mathbf{S}_i \cdot \mathbf{S}_{i+1} = J \sum_{i=1}^{N-1} (S_i^x S_{i+1}^x + S_i^y S_{i+1}^y + S_i^z S_{i+1}^z).
$$ 

In the DSSF, the state $|\psi_0\rangle$ is the ground state of the AFH model. The DSSF is directly related to what can be measured in inelastic neutron scattering experiments, and we want to use MPS methods to provide the theory predicitions corresponding to the experimental results shown in {numref}`fig:lake2013_2`.

```{figure} ../week1/images/AFH_experiment.png
---
name: fig:lake2013_2
width: 80%
align: center
---
Data from B. Lake *et al* {cite}`Lake2013`. The left panel shows the experimental data for the dynamical structure factor of $\text{KCuF}_3$. The right panel shows the theoretical prediction for the Heisenberg model obtained using Bethe Ansatz.
```


## Simplifying the DSSF

Before setting up the calculation using MPS, we are going to simplify the DSSF expression in {eq}`eq:dsf_2`. Firstly, we are only interested in the diagonal components $S^{\alpha\alpha}$. Furthermore, due to the symmetry of the Hamiltonian, we can consider only $S^{zz}$. Then we can use an approximation that the system is translation invariant. We ideally would like to compute the DSSF for an infinite translation invariant system, but instead our MPS code will use a (large) finite system size $N$. This approximation gets better as we increase $N$. With this approximation, we can write the DSSF as 

$$
\begin{aligned}
S^{\alpha\alpha}(q, \omega) &= \frac{1}{2\pi} e^{-iqN/2}\sum_{r} e^{iqr} \int_{-\infty}^{\infty} e^{i\omega t} \langle \psi_0 | S^\alpha_r(t) S^\alpha_{N/2}(0) | \psi_0 \rangle \;dt,\\
&= \frac{1}{2\pi} \sum_{r=-N/2}^{N/2-1} e^{iqr} \int_{-\infty}^{\infty} e^{i\omega t} \langle \psi_0 | S^\alpha_{r+N/2}(t) S^\alpha_{N/2}(0) | \psi_0 \rangle \;dt,
\end{aligned}
$$ (eq:dsf_translation_invariant)

where we have fixed the second spin at the centre of the system to minimise the finite size effects. Next we don't want to have to compute the integral to negative times. We can rewrite the integral to get only one side of the limit. We can rewrite the integral as

$$
\begin{aligned}
&\int_{0}^{\infty} e^{i\omega t} \langle \psi_0 | S^\alpha_{r-N/2}(t) S^\alpha_{N/2}(0) | \psi_0 \rangle + \int_{-\infty}^{0} e^{i\omega t}\langle \psi_0 | S^\alpha_{r+N/2}(t) S^\alpha_{N/2}(0) | \psi_0 \rangle \;dt \\
&= \int_{0}^{\infty} e^{i\omega t} \langle \psi_0 | S^\alpha_{r-N/2}(t) S^\alpha_{N/2}(0) | \psi_0 \rangle + \int_{0}^{\infty} e^{-i\omega t}\langle \psi_0 | S^\alpha_{r+N/2}(-t) S^\alpha_{N/2}(0) | \psi_0 \rangle \;dt \\
&= \int_{0}^{\infty} e^{i\omega t} \langle \psi_0 | S^\alpha_{r-N/2}(t) S^\alpha_{N/2}(0) | \psi_0 \rangle + \int_{0}^{\infty} \left(e^{i\omega t}\langle \psi_0 | S^\alpha_{N/2}(t) S^\alpha_{r+N/2}(0) | \psi_0 \rangle \right)^* \;dt \\
&= \int_{0}^{\infty} e^{i\omega t} \langle \psi_0 | S^\alpha_{r-N/2}(t) S^\alpha_{N/2}(0) | \psi_0 \rangle + \int_{0}^{\infty} \left(e^{i\omega t}\langle \psi_0 | S^\alpha_{r+N/2}(-t) S^\alpha_{N/2}(0)| \psi_0 \rangle \right)^* \;dt \\
&=\int_{0}^{\infty} 2\text{Re}\left[e^{i\omega t} \langle \psi_0 | S^\alpha_{r+N/2}(t) S^\alpha_{N/2}(0) | \psi_0 \rangle\right] \;dt.
\end{aligned}
$$

In the second line, we flipped the integral limits, and sent $t\rightarrow -t$. In the third line, we used the fact that $\langle S^\alpha_{j}(-t) S^\alpha_{k}(0) \rangle = \langle S^\alpha_{k}(0) S^\alpha_{j}(-t) \rangle^* = \langle S^\alpha_{k}(t) S^\alpha_{j}(0) \rangle^*$. Then in the fourth line we again use the translation invariance to shift the indices, along with the fact that we sum over all $r$ and so can flip the sign of $r$. The final line simply follows from $2\text{Re}(z) = z + z^*$. We therefore get the integral

$$
S^{\alpha\alpha}(q, \omega) = \frac{1}{2\pi} \sum_{r} e^{iqr} \int_{0}^{\infty} 2\text{Re}\left[e^{i\omega t} \langle \psi_0 | S^\alpha_{r+N/2}(t) S^\alpha_{N/2}(0) | \psi_0 \rangle\right] \;dt.
$$

Finally, we replace the integral by a sum over discrete time steps $\Delta t$ from our simulation to get

$$
S^{\alpha\alpha}(q, \omega) \approx 2\sum_{r=-N/2}^{N/2-1} \sum_{m=0}^{M} e^{iqr} \text{Re}\left[e^{i\omega m \Delta t} \langle \psi_0 | S^\alpha_{r+N/2}(t) S^\alpha_{N/2}(0) | \psi_0 \rangle\right]
$$ (eq:dsf_final)

### Removing Gibbs oscillations

```{figure} images/gibbs.jpg
---
name: fig:gibbs
width: 85%
align: center
---

Removing Gibbs oscillations using a window function. Left shows the window function, a correlator, and the windowed correlator. The right shows the Fourier transform of the correlator and windowed correlator. 
``` 


In replacing the integral over time with a finite sum over discrete time steps, we have introduced two sources of error. The first is due to the finite time step $\Delta t$, and the second is due to the finite time window $t = M\Delta t$. The finite time step error can be reduced by decreasing $\Delta t$. The finite time window introduces a sharp cutoff to the integral. This sharp cutoff can introduce Gibbs oscillations, as can be seen in {numref}`fig:gibbs`. To remove these oscillations, we can use a window function to introduce a smooth cutoff to the integral. There are many options for window functions. We will use a simple cosine squared window function (shown in {numref}`fig:gibbs`)

$$
g(m) = \cos^2\left(\frac{\pi m}{2M}\right),
$$ (eq:window)

where the argument is chosen such that it is zero at $t=M\Delta t$. As we can see in {numref}`fig:gibbs`, the finite time window introduces oscillations in the Fourier transform of the correlator. By taking the Fourier transform of the windowed correlator, we can see that the oscillations are removed and the signal is smoothed.

We can now write the final expression for the DSSF as

$$
S^{\alpha\alpha}(q, \omega) \approx 2\sum_{r=-N/2}^{N/2-1} \sum_{m=0}^{M} e^{iqr} \text{Re}\left[e^{i\omega m \Delta t} \langle \psi_0 | S^\alpha_{r+N/2}(t) S^\alpha_{N/2}(0) | \psi_0 \rangle \right] g(m)
$$ 

where we have introduced the window function $g(m)$ to smooth the correlator.

```{note}
It is a bit unfortunate that the window function weights our correlator such that the longest times have the least impact on the final result. However, the later times are those that are hardest to compute, so it feels like we are throwing away that hard work. In practice, it is common to use a simple extrapolation method to extend the results of the simulation to longer times. These extropolated results cannot be fully trusted, but since they are strongly suppressed by the window function, they do not have a large impact on the final result. The time steps that we put hard work into calculating will now have a larger impact on the final result.

For simplicity, we won't introduce any extrapolation methods in this course.
```

## Setting up the calculation

Now that we have figured out how to correctly compute the Fourier transform, all that remains is to use MPS to compute the correlation functions

$$
C^{zz}_j(t) = \langle \psi_0 | S^z_{j}(t) S^z_{N/2}(0) | \psi_0 \rangle.
$$

Again, it is worth the effort to manipulate this expression to simplify our calculation. First let us reveal the time evolution of the operators

$$
\begin{aligned}
C^{zz}_j(t) &= \langle \psi_0 | e^{iHt} S^z_{j} e^{-iHt} S^z_{N/2} | \psi_0 \rangle, \\
&= e^{iE_0 t} \langle \psi_0 | S^z_{j} e^{-iHt} S^z_{N/2} | \psi_0 \rangle, \\
&= e^{iE_0 t} \left(\langle \psi_0 | S^z_{j} e^{-iHt/2}\right) \left( e^{-iHt/2} S^z_{N/2} | \psi_0 \rangle \right), \\
&= 4e^{iE_0 t} \left(\langle \psi_0 | \sigma^z_{j} e^{-iHt/2}\right) \left( e^{-iHt/2} \sigma^z_{N/2} | \psi_0 \rangle \right).
\end{aligned}
$$

In the second line, we have used the fact that $|\psi_0\rangle$ is the ground state of $H$, i.e. it is an eigenstate with eigenvalue $E_0$. In the third line, we have split the time evolution into two parts. This allows us to write the correlator as the overlap between two states that have been evolved only to time $t/2$. This is useful because computing the time evolution for long times is costly due to the growth of entanglement, and so we can save some computational effort. Finally, we rewrite in terms of the Pauli operators $\sigma^z$, instead of the spin operators $S^z = \sigma^z$. This is because the Pauli operators are unitary, and so we can act with them on the MPS without changing the canonical form.




````{admonition} Algorithm: Computing the correlator
Let us now wrtie down the steps to compute the correlator $C^{zz}_j(t)$.

1. Compute the ground state $|\psi_0\rangle$ of the AFH model using DMRG.

2. Compute the state $|\psi_{N/2}\rangle = \sigma^Z_{N/2} |\psi_0\rangle$. and the state $|\psi_j\rangle = \sigma^Z_{j} |\psi_0\rangle$.

3. Evolve the states using TEBD to time $t/2$. More explicitly, we evolve $|\psi_{N/2}(t/2)\rangle = e^{-iHt/2} |\psi_{N/2}\rangle$ and $|\psi_j(t/2)\rangle = e^{iHt/2} |\psi_j\rangle$. Note that evolve "backwards" in time for $|\psi_{j}(t/2)\rangle$, which is done by setting $\Delta t \rightarrow -\Delta t$ in the TEBD algorithm.

4. Compute the overlap between $|\psi_{N/2}(t/2)\rangle$ and $|\psi_j(t/2)\rangle$ to get the correlator $C^{zz}_j(t)$. More precisely, compute $C^{zz}_j(t) = 4 e^{iE_0t}\langle \psi_{j}(t/2) | \psi_{N/2}(t/2) \rangle$.

Note that when doing the time evolution we can evolve one of the two states at a time. That is, evolve $|\psi_{N/2}(t/2)\rangle$ forward one step, and then compute the overlap. Then evolve $|\psi_j(t/2)\rangle$ backwards one step, and then compute the overlap. This way we only need to do the evolution $M$ times and not $2M$ times.

````

To compute the correlator, there are two operations that we need to perform that we have not yet implemented: applying an operator to a state, and computing the overlap between two states. These are simple operations on the MPS that we will add to our code.

### Applying an operator to a state


```{figure} images/apply_local.jpeg
---
name: fig:apply_local
width: 45%
align: center
---

Applying a local unitary to an MPS can be done with a single local contraction. The operator is applied to the physical index of the MPS.
``` 

The first new operation is acting on an MPS with a local operator. In our case we only need to act with the $\sigma^z$ operator. This is local and unitary. Because the operator is unitary, it does not change the canonical form, meaning we can do this operation with a single local contraction without needing to move the centre.


````{admonition} Code: Applying an operator to a state

```python
## file: src/mps.py


```
```{figure} images/apply_contraction.jpeg
---
name: fig:apply_local
width: 60%
align: center
---

Indexing for the local contraction of an operator applied to an MPS tensor.
``` 


````


### Computing the overlap

```{figure} images/overlap.jpeg
---
name: fig:overlap
width: 75%
align: center
---

Tensor network diagram for the overlap between two MPS.
``` 


The second new operation is computing the overlap between two states. Since the states we computing the overlap between are different states, we cannot take advantage of the canonical form, and therefore there is also no need to move the centre. 


---

## References

```{bibliography}
:filter: docname in docnames
```