# Bringing It All Together

At this point we have learnt all the aspects of MPS necessary to tackle the problem we set out at the start: computing the dynamical spin structure factor for the antiferromagnetic Heisenberg model. We have learnt about the MPS representation of quantum states, the DMRG algorithm for finding ground states, and the TEBD for simulating time evolution. It is now time to bring all these elements together.

## Dynamical Spin Structure Factor

Let us restate the problem that we are tackling with our MPS code. We want to compute the *dynamical spin structure factor* (DSSF)

$$
S^{\alpha\beta}(q, \omega) = \frac{1}{2\pi N}\sum_{j,k} e^{iq(j-k)} \int_{-\infty}^{\infty} e^{i\omega t} \langle \psi_0 | S^\alpha_j(t) S^\beta_k(0) | \psi_0 \rangle \;dt,
$$ (eq:dsf)

for the antiferromagnetic Heisenberg (AFH) model described by the Hamiltonian

$$
H = J \sum_{i=1}^{N-1} \mathbf{S}_i \cdot \mathbf{S}_{i+1} = J \sum_{i=1}^{N-1} (S_i^x S_{i+1}^x + S_i^y S_{i+1}^y + S_i^z S_{i+1}^z).
$$ (eq:heisenberg)

In the DSSF, the state $|\psi_0\rangle$ is the ground state of the AFH model. The DSSF is directly related to what can be measured in inelastic neutron scattering experiments, and we want to use MPS methods to provide the theory predicitions corresponding to the experimental results shown in {numref}`fig:lake2013`.

```{figure} ../week1/images/AFH_experiment.png
---
name: fig:lake2013
width: 80%
align: center
---
Data from B. Lake *et al* {cite}`Lake2013`. The left panel shows the experimental data for the dynamical structure factor of $\text{KCuF}_3$. The right panel shows the theoretical prediction for the Heisenberg model obtained using Bethe Ansatz.
```


## Simplifying the DSSF

Before setting up the calculation using MPS, we are going to simplify the DSSF expression in {eq}`eq:dsf`. Firstly, we are only interested in the diagonal components $S^{\alpha\alpha}$. Furthermore, due to the symmetry of the Hamiltonian, we can consider only $S^{zz}$. Then we can use an approximation that the system is translation invariant. We ideally would like to compute the DSSF for an infinite translation invariant system, but instead our MPS code will use a (large) finite system size $N$. This approximation gets better as we increase $N$. With this approximation, we can write the DSSF as 

$$
\begin{aligned}
S^{\alpha\alpha}(q, \omega) &= \frac{1}{2\pi} e^{-iqN/2}\sum_{r} e^{iqr} \int_{-\infty}^{\infty} e^{i\omega t} \langle \psi_0 | S^\alpha_r(t) S^\alpha_{N/2}(0) | \psi_0 \rangle \;dt,\\
&= \frac{1}{2\pi} \sum_{r=-N/2}^{N/2-1} e^{iqr} \int_{-\infty}^{\infty} e^{i\omega t} \langle \psi_0 | S^\alpha_{r+N/2}(t) S^\alpha_{N/2}(0) | \psi_0 \rangle \;dt,
\end{aligned}
$$ (eq:dsf_translation_invariant)

where we have fixed the second spin at the centre of the system to minimise the finite size effects. Next we don't want to have to compute the integral to negative times. We can rewrite the integral to get only one side of the limit. We can rewrite the integral as

$$
\begin{aligned}
&= \int_{0}^{\infty} e^{i\omega t} \langle \psi_0 | S^\alpha_{r-N/2}(t) S^\alpha_{N/2}(0) | \psi_0 \rangle + \int_{-\infty}^{0} e^{i\omega t}\langle \psi_0 | S^\alpha_{r+N/2}(t) S^\alpha_{N/2}(0) | \psi_0 \rangle \;dt \\
&= \int_{0}^{\infty} e^{i\omega t} \langle \psi_0 | S^\alpha_{r-N/2}(t) S^\alpha_{N/2}(0) | \psi_0 \rangle + \int_{0}^{\infty} e^{-i\omega t}\langle \psi_0 | S^\alpha_{r+N/2}(-t) S^\alpha_{N/2}(0) | \psi_0 \rangle \;dt \\
&= \int_{0}^{\infty} e^{i\omega t} \langle \psi_0 | S^\alpha_{r-N/2}(t) S^\alpha_{N/2}(0) | \psi_0 \rangle + \int_{0}^{\infty} \left(e^{i\omega t}\langle \psi_0 | S^\alpha_{N/2}(t) S^\alpha_{r+N/2}(0) | \psi_0 \rangle \right)^* \;dt \\
&= \int_{0}^{\infty} e^{i\omega t} \langle \psi_0 | S^\alpha_{r-N/2}(t) S^\alpha_{N/2}(0) | \psi_0 \rangle + \int_{0}^{\infty} \left(e^{i\omega t}\langle \psi_0 | S^\alpha_{r+N/2}(-t) S^\alpha_{N/2}(0)| \psi_0 \rangle \right)^* \;dt \\
&=\int_{0}^{\infty} 2\text{Re}\left[e^{i\omega t} \langle \psi_0 | S^\alpha_{r+N/2}(t) S^\alpha_{N/2}(0) | \psi_0 \rangle\right] \;dt.
\end{aligned}
$$

In the second line, we flipped the integral limits, and sent $t\rightarrow -t$. In the third line, we used the fact that $\langle S^\alpha_{j}(-t) S^\alpha_{k}(0) \rangle = \langle S^\alpha_{k}(0) S^\alpha_{j}(-t) \rangle^* = \langle S^\alpha_{k}(t) S^\alpha_{j}(0) \rangle^*$. Then in the fourth line we again use the translation invariance to shift the indices, along with the fact that we sum over all $r$ and so can flip the sign of $r$. The final line simply follows from $2\text{Re}(z) = z + z^*$. We therefore get the integral

$$
S^{\alpha\alpha}(q, \omega) = \frac{1}{2\pi} \sum_{r} e^{iqr} \int_{0}^{\infty} 2\text{Re}\left[e^{i\omega t} \langle \psi_0 | S^\alpha_{r+N/2}(t) S^\alpha_{N/2}(0) | \psi_0 \rangle\right] \;dt.
$$ (eq:dsf_final)

Finally, we replace the integral by a sum over discrete time steps $\Delta t$ from our simulation to get

$$
S^{\alpha\alpha}(q, \omega) \approx 2\sum_{r=-N/2}^{N/2-1} \sum_{m=0}^{M} e^{iqr} \text{Re}\left[e^{i\omega m \Delta t} \langle \psi_0 | S^\alpha_{r+N/2}(t) S^\alpha_{N/2}(0) | \psi_0 \rangle\right]
$$ (eq:dsf_final)

### Removing Gibbs oscillations

```{figure} images/gibbs.jpg
---
name: fig:gibbs
width: 80%
align: center
---

Removing Gibbs oscillations using a window function. Left shows the window function, a correlator, and the windowed correlator. The right shows the Fourier transform of the correlator and windowed correlator. 
``` 


In replacing the integral over time with a finite sum over discrete time steps, we have introduced two sources of error. The first is due to the finite time step $\Delta t$, and the second is due to the finite time window $t = M\Delta t$. The finite time step error can be reduced by decreasing $\Delta t$. The finite time window effective introduces a sharp cutoff to the integral. This sharp cutoff can introduce Gibbs oscillations, as can be seen in {numref}`fig:gibbs`. To remove these oscillations, we can use a window function to introduce a smooth cutoff to the integral. There are many options for window functions. We will use a simple cosine squared window function (shown in {numref}`fig:gibbs`)

$$
g(m) = \cos^2\left(\frac{\pi m}{2M}\right).
$$ (eq:window)

As we can see in {numref}`fig:gibbs`, the finite time window introduces oscillations in the Fourier transform of the correlator. By taking the Fourier transform of the windowed correlator, we can see that the oscillations are removed and the signal is smoothed.

```{note}
It is a bit unfortunate that the window function weights our correlator such that the longest times have the least impact on the final result. However, the later times are those that are hardest to compute, so it feels like we are throwing away that hard work. In practice, it is common to use a simple extrapolation method to extend the results of the simulation to longer times. These extropolated results cannot be fully trusted, but since they are strongly suppressed by the window function, they do not have a large impact on the final result. The time steps that we put hard work into calculating will now have a larger impact on the final result.

For simplicity, we won't introduce any extrapolation methods in this course.
```