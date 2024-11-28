# Density Matrix Renormalization Group (DMRG)

Now that we have introduced MPOs, we are in a position to discuss the Density Matrix Renormalization Group (DMRG) algorithm. DMRG is a numerical method for finding the ground state of a 1D quantum system. It was first introduced by Steven White in 1992 {cite}`White1992` and has since become one of the most powerful methods for studying 1D quantum systems. The algorithm is iterative and involves sweeping through the system, updating the tensors in the MPS representation. We will explain a variant of the DMRG algorithm that uses a two-site update. This is the simplest variant of DMRG since it naturally allows the bond dimension to grow as we sweep through the system. There is also a one-site update variant, which we will not discuss in this course.

The basis of the DMRG algorithm is to locally update the tensors in the MPS representation in order to minimise $\langle \psi | H |\psi\rangle$. This is done by considering two sites at a time and finding the ground state for an effective Hamiltonian on these two sites. The process is summarised in {numref}`dmrg_blocks`. Let us consider updating the tensors on site $n$ and $n+1$. We can split the computation of the energy in to several parts. First, after moving the centre to site $n$, we can contract all the tensors to the left of site $n$, to construct a rank-3 tensor $L^{[n-1]}$, which we refer to as the left environment. Similarly, we can construct the right environment $R^{[n+2]}$ by contracting all tensors to the right of site $n+1$. In between the left and right environments, we have the tensors from the MPO, the MPS and its conjugate on site $n$ and $n+1$. We can further combine the tensors from the MPS into $\Theta$, as shown in {numref}`dmrg_blocks`. The goal is then to update the tensor $\Theta$ in order to minimise the value from contracting the entire diagram. 



```{figure} images/dmrg_blocks.jpeg
---
name: fig:dmrg_blocks
width: 100%
align: center
---

Blocks of tensor contractions to set up the DMRG algorithm with a two-site update. The diagram computes the energy of the MPS $\langle \psi | H | \psi \rangle$. When updating sites $n$ and $n+1$, we can contract all the tensors to the left of $n$ into a left environment, and similarly can construct a right environment from all those right of $n+1$. The DMRG algorithm aims to locally change the grouped tensor $\Theta$ in order to minimise the energy.
```

In order to update the tensor $\Theta$, we construct an effective Hamiltonian for the sites $n$ and $n+1$. This is done by simply removing $\Theta$ and its conjugate from the diagram in {numref}`dmrg_blocks` to get the rank-8 tensor shown in {numref}`h_effective`. We can then reshape this tensor into a matrix. The problem of minimising the value of the entire tensor network contraction is equivalent to finding a new $\tilde{\Theta}$, which is the ground state of $H_\text{eff}$.



```{figure} images/h_effective.jpeg
---
name: fig:h_effective
width: 70%
align: center
---

The energy minimisation problem when using a two-site update is equivalent to finding the ground state of the effective Hamiltonian shown in this diagram. The resulting rank-8 tensor can be reshaped into a matrix to use standard eigensolvers.
```

```{note}
Solving the local eigenproblem is the most costly part of the DMRG algorithm. In this course we will solve it using a sparse iterative solver `scipy.sparse.linalg.eigsh`. This allows us to find only the lowest eigenvalue and eigenstate without solving the full eigenproblem.

I practice, for performance, the eigensolver should be replaced by a custom iterative solver. You can then approximately solve this eigenproblem with just a few iterations, using the previous solution as a starting point. As few as 2 iterations can be sufficient for the DMRG algorithm to converge since you are sweeping through the system and each tensor will be updated multiple times. Reducing the time spent on each eigenproblem dramatically increases the performance. We won't do this in this course.
```


After finding the ground state of the effective Hamiltonian, we want to split the new $\tilde{\Theta}$ tensor into the tensors of our MPS. As always, we do this by performing a truncated SVD, as shown in {numref}`fig:svd`. First we reshape the vector into a matrix, with the left virtual and physical legs grouped, and similarly for the right. We can then perform the truncated SVD and combine the Schmidt values into the right tensor. For a left to right sweep, we will start with the centre on site $n$ and move it to $n+1$ after the SVD. For the right to left sweep we will do the opposite and start with the centre on site $n+1$ and move it to $n$. By doing this we never have to explicitly move the centre, we make the move as a biproduct of having the do SVD. All we will have to do is move the centre at the very start of the DRMG algorithm.


```{figure} images/svd.jpeg
---
name: fig:svd
width: 80%
align: center
---

We can extract the MPS tensors from the new $\tilde{\Theta}$ tensor by performing truncated SVD.
```

The DMRG algorithm then proceeds by sweeping through the system from left to right and back again, updating two sites at a time until it converges.





```{figure} images/update_L.jpeg
---
name: fig:update_L
width: 45%
align: center
---

???
```


```{figure} images/h_contraction.jpeg
---
name: fig:h_contraction
width: 75%
align: center
---

???
```