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


