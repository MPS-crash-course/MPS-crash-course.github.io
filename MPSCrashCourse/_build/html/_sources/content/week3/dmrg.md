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

The DMRG algorithm then proceeds by sweeping through the system from left to right and back again, updating two sites at a time until it converges. In the process we will have to compute new left and right environments. Instead of computing the contract of the entire left/right environment every time, the environment can be updated iteratively, as shown in {numref}`fig:update_L`. We can use the `get_slice` method of the MPO that we defined last time to compute the slice using the updated MPS. When sweeping from left to right, we can contract this slice with $L^{[n-1]}$ we get $L^{[n]}$. This is because none of the tensors in the MPS to the left of site $n$ have changed. We can update the right environment similarly when sweeping from right to left. Finally, we can store a list of all the environments to the left and right of the current site, so that we can reuse them when we sweep back. More explicitly, as we sweep left to right, we add the new $L^{[n]}$ to the list of left environments, and then we remove $R^{[n+2]}$ from the list of right environments since it involves tensors that will be updated in the next step.



```{figure} images/update_L.jpeg
---
name: fig:update_L
width: 45%
align: center
---

Computing the left environement tensor $L^{[n]}$ from $L^{[n-1]}$ by contracting with the next slice. This avoids computing the full contraction of the left environment each time.
``` 

````{admonition} Algorithm: DMRG

Let us summarise the steps of the DMRG algorithm with a two-site update:

1. Start by creating a product state MPS and move the centre to the far right of the MPS.

2. Iterative construct the list of left environments (leaving the last two sites of the chain). Start with $L^{[0]} = 1$ and then for each site $n$ contract the slice of the MPO with $L^{[n-1]}$ to get $L^{[n]}$.

3. Contruct the list of right environements. This will contain a single element $R^{[N+1]} = 1$.

4. Sweep from right to left (repeat until we end at the first two sites):
    1. Construct effective Hamiltonian for sites $n$ and $n+1$.
    2. Find the ground state of the effective Hamiltonian.
    3. Perform truncated SVD to update the MPS tensors (moving the centre one step to the left).
    4. Update the right environment tensor $R^{[n+1]}$ by contracting the slice of the MPO with $R^{[n+2]}$ and add to list of right environments.
    5. Remove $L^{[n-1]}$ from the list of left environments.

5. Sweep from left to right (repeat until we end at the last two sites):
    1. Construct effective Hamiltonian for sites $n$ and $n+1$.
    2. Find the ground state of the effective Hamiltonian.
    3. Perform truncated SVD to update the MPS tensors (moving the centre one step to the right).
    4. Update the left environment tensor $L^{[n]}$ by contracting the slice of the MPO with $L^{[n-1]}$ and add to list of left environments.
    5. Remove $R^{[n+2]}$ from the list of right environments.

6. Measure and store the energy with respect to the current MPS. Repeat steps 4 and 5 until we hit the maximimum number of sweeps.

It would be best to add a convergence criteria to stop the algorithm earlier if it converges, but we won't do that in this course.

````



````{admonition} Code: DMRG Algorithm

```python
## file: src/dmrg.y

import numpy as np
from .mpo import *
from .mps import *
from .svd import svd_truncated

import scipy.sparse as sp


def dmrg(psi, H_mpo, chiMax, tol=1E-12, nSweeps=5):
    """
    Perform density matrix renormalization group (DMRG) for a matrix product state (MPS).

    Parameters
    ----------
    psi : MPS
        Initial matrix product state.
    H_mpo : MPO
        Hamiltonian as a matrix product operator.
    chiMax : int
        Maximum bond dimension.
    tol : float, optional
        Convergence tolerance.
    nSweeps : int, optional
        Number of DMRG sweeps.

    Returns
    -------
    psi : MPS
        Ground state matrix product state.
    E : float
        Ground state energy.
    """
    
    L = psi.L
    assert L == H_mpo.L, "MPS and MPO sizes do not match."
    
    # move centre to right end (won't do anything for product state)
    psi.move_centre_to(L-1)  

    # build left environment list
    L_envs = [np.array([1]).reshape((1,1,1))]
    for i in range(L-2):
        ## YOUR CODE HERE ##

    # build right environment
    R_envs = [np.array([1]).reshape((1,1,1))]

    E_list = []
    for _ in range(nSweeps):

        # sweep right to left in dmrg
        for i in range(L-2, 0, -1):
            psi.move_centre_to(i+1)  # Shouldn't actually do anything (but just to be safe)

            H_block, chi_left, chi_right = construct_H_block(H_mpo, i, L_envs, R_envs)

            ## YOUR CODE HERE ##
        
        # sweep left to right in dmrg
        for i in range(L-2):
            psi.move_centre_to(i)

            H_block, chi_left, chi_right = construct_H_block(H_mpo, i, L_envs, R_envs)

            ## YOUR CODE HERE ##

        # compute final energy
        E_list.append(H_mpo.expectation(psi))

    return psi, E_list


def construct_H_block(H_mpo, i, L_envs, R_envs):
    ## YOUR CODE HERE ##

    chi_left, d2, d3, chi_right, d5, d6, d7, d8 = H_block.shape
    H_block = H_block.reshape(chi_left*d2*d3*chi_right, d5*d6*d7*d8)
    return H_block, chi_left, chi_right
```

In order to contrust the $H_\text{eff}$ Hamiltonian we have defined a function called `construct_H_block`. This should perform the contraction shown in {numref}`fig:h_contraction`.

```{figure} images/h_contraction.jpeg
---
name: fig:h_contraction
width: 75%
align: center
---

Contractions to compute the $H_\text{eff}$ tensor. The red indices correspond to the first input to `tensordot`, and green are the second. The final step is a transpose of the indices.
```


````






## Testing the DMRG algorithm

The DMRG algorithm is significantly more complicated than the TEBD algorithm, so it is important that we test that it is working as expected. To do this we should compare the energy to results from exact diagonalization for small system sizes. You should write code to plot the absolute error as a function of DMRG sweeps (iterations), for different bond dimensions, to get the results shown in {numref}`fig:dmrg_test`. 


```{figure} images/dmrg_test.png
---
name: fig:dmrg_test
width: 60%
align: center
---

Testing the DMRG algorithm on the Antiferromagnetic Heisenberg model against exact diagonalization for different maximum bond dimensions. The absolute error in the energy is shown. Results are for system size $L=10$, tolerance = $10^{-14}$.
```

We can see that the DMRG algorithm very quickly converges in energy up to a limit in accuracy set by the maximum bond dimension. For this system size ($L=10$), the algorithm has converged after only 3 full sweeps. By increasing the bond dimension we increase the accuracy of the final result. By bond dimension 32 ($=2^5$), we are able to exactly represent the state for $L=10$ states, and can get the correct energy up to machine precision.

```{note}
Now we have set up the DMRG algorithm and tested it, we can find the ground state for larger systems, beyond the reach of exact diagonalization. When doing this, we need to check that the energy is converged in the bond dimension. 

The DMRG algorithm is particularly powerful for finding ground states of *gapped* one-dimensional systems, since the ground states satisfy an *area-law*. This means that as we increase the system size, we don't need to increase the bond dimension. However, the Hamiltonian we have chosen to look at in this course is actually *gapless*, and so we will need to increase the bond dimension with system size to achieve a comparable accuracy. This makes it one of the harder cases for DMRG. Thankfully, we will be able to access a system size of 100 with only a moderate bond dimension of 16 for our final simulation next week and get sufficiently accurate results.

```


````{admonition} Exercise: Large system entanglement

It would be good to get some practice with using the DMRG and to use it on system size beyond ED, so let us write an exercise. I added the file `exercises/dmrg_exercise.py`. In this exercise, find the ground state of the Antiferromagnetic Heisenberg model for system size $L=100$ and bond dimension $\chi=16$ and 5 sweeps. From this state, compute the entanglement entropy with respect to every bond in the chain, and plot entropy vs bond. Repeat this with system size $L=101$. You should find the results shown in {numref}`fig:entanglement_profile`.

```{figure} images/entanglement_profile.png
---
name: fig:entanglement_profile
width: 90%
align: center
---

The entanglement entropy profile for the ground state of the Antiferromagnetic Heisenberg chain for system size $L=100$ and $L=101$, computed using DMRG.
```

From these results we can see behaviour related to the gaplessness of the model. First of all, there is a characteristic "arch" shape to the entanglement profile. This contrasts area-law, which would give a flat profile, and volume-law, which would give a triangle shaped profile. Instead, the centre of the profile will increase logarithmically with system size. The second feature is the dramatic difference in the profile when changing the system size by just one site. Gapless systems do not have a finite correlation length and instead have polynomially decaying correlations. This is revealed here by a sensitivity to the boundary of the system. 

````
