# Entanglement and Truncation

One of the central concepts in MPS is entanglement. Entanglement is a measure of how correlated two parts of a quantum system are. The entanglement spectrum (Schmidt values) also allows us to controllably truncate the MPS representation.

## The Schmidt decomposition

```{figure} images/schmidt.jpeg
---
name: fig:schmidt
width: 50%
align: center
---

The Schmidt decomposition of a quantum state with respect to a bipartition of the system into two contiguous parts $A$ and $B$, written in tensor network diagrams. The Schmidt values are contained in a diagonal matrix $\Lambda$, and can be used to computed the entanglement entropy between the two parts.
```

Given a pure quantum state $|\psi\rangle$, and a bipartition of the system into two parts $A$ and $B$, the Schmidt decomposition of $|\psi\rangle$ is given by

$$
|\psi\rangle = \sum_i \lambda_i |i_A\rangle \otimes |i_B\rangle,
$$

where $\{|i_A\rangle\}$ and $\{|i_B\rangle\}$ are orthonormal bases for parts $A$ and $B$, respectively, and $\lambda_i$ are the Schmidt values. The Schmidt values are non-negative and are related to the entanglement of the state. In particular, the entanglement between the two subsystems can be quantified by the von Neumann entanglement entropy

$$
S = -\sum_i \lambda_i^2 \log \lambda_i^2.
$$

The MPS form a quantum state naturally gives us access to the Schmidt values. In particular, the Schmidt decomposition of a state can be written as the tensor network diagram shown in {numref}`fig:schmidt`. However, because we have chosen not to introduce the central tensors, the Schmidt values are not directly accessible. We therefore need to perform SVD in order to extract the Schmidt values, as shown in {numref}`fig:extract_schmidt`. If we want to measure the entanglement entropy on the bond between sites $n$ and $n+1$, then we need to first move the centre to site $n$ and then perform the SVD.

```{figure} images/extract_schmidt.jpeg
---
name: fig:extract_schmidt
width: 80%
align: center
---

In our code we don't explicitly expose the Schmidt values and so need to extract them from the central tensor. This is done by contracting with a neighbouring tensor and then performing SVD. Here we have contracted to the right, but could have equivalently contracted to the left, assuming that all tensors are in left (right) canonical form to the left (right) of the centre tensor.
```

```{note}

Note that the singular values from SVD are not necessarily the Schmidt values. For instance, if we did not move the centre and instead combined two sites into one then performed the SVD, the singular values would not be the Schmidt values. In fact, we would have no physical interpretation of the singular values in this case. This is because the states to the left and right of this bond are not orthonormal. Orthonormality is only guaranteed when the centre is at the site of interest.

For us a bipartition is always defined with respect to a bond. The subsystems $A$ and $B$ are the contiguous sites to the left and right of the bond, respectively. While other bipartitions are possible, it is not easy to extract the Schmidt values from the MPS representation in those cases.

```


````{admonition} Code: Entanglement entropy

It will be useful for us to measure the bipartite entanglement entropy of our state, so let's add this to our MPS class.

```python
## file: src/mps.py

class MPS:

    ## PREVIOUS CODE OMITTED ##

    def entropy(self, site):
        """
        Compute the bipartite entanglement entropy of the bond between site and site+1.
        """
        # Move the centre to the site of interest
        self.move_centre_to(site)

        ## YOUR CODE HERE ##

        return -np.sum(S**2 * np.log(S**2))


```

````


## Truncation

```{figure} images/truncated_svd.jpeg
---
name: fig:truncated_svd
width: 80%
align: center
---

We will encounter tensors on two sites in the TEBD and DMRG algorithms. In order to split this tensor into two tensors, we need to perform SVD. We can truncate this SVD to limit the bond dimension of our MPS by keeping only the largest Schmidt values and the corresponding columns of $U$ and rows of $V^\dagger$.
```

The Schmidt values provide us with a controlled way to truncate the MPS representation. The idea is to keep only the largest Schmidt values and discard the rest. This is equivalent to keeping only the most important correlations in the state. As shown in {numref}`fig:truncated_svd`, in the TEBD and DMRG algorithms, we will end up with a tensor on two sites. We then want to perform SVD to split this tensor into two tensors, one on each site. In order to limit the size of the tensors in our MPS, we can truncate this SVD. This is done by keeping only the largest $\chi$ singular values, only the $\chi$ columns of $U$ and the $\chi$ rows of $V^\dagger$.

To determine how many singular values to keep, we can introduce two thresholds. The first is a threshold on the accuracy of the truncation, $\epsilon$. We keep all singular values such that

$$
\sum_{i=1}^{\chi} \lambda_i^2 \geq 1 - \epsilon.
$$

Because the Schmidt values are related to the Schmidt decomposition of the state (so long as this tensor is at the centre), the infidelity between the correct state and the truncated state is given by $\epsilon$. This is one of the powerful aspects of MPS methods, that it is a controlled approximation, and we can keep track of the errors accumulated in our simulations. The second threshold is to set a maximum allowable $\chi = \chi_\text{max}$. This is useful when we want to limit the computational cost of the simulation since the required $\chi$ to satisfy our accuracy threshold can be very large, particularly for highly entangled states. When we truncate the Schmidt values, the state will no longer be normalized since $\sum_{i}^{\chi} \lambda_i^2 < 1$. We therefore also need to renormalize the Schmidt values.


````{admonition} Code: Truncated SVD

Here I provide code for the truncated SVD. This function takes a matrix $M$, so you may need to reshape your tensor before and after applying this SVD. The function returns the truncated $U$, $S$ and $V^\dagger$ matrices, where the bond dimension is determined by the accuracy threshold $\epsilon$, and the maximum bond dimension $\chi_\text{max}$. We will add this code to a new file called `svd.py` in the `src` folder.

```python
## file: src/svd.py

import numpy as np
import numpy.linalg as la

def svd_truncated(M, chiMax, threshold):

    U, S, Vdg = la.svd(M, full_matrices=False)
    
    if (chiMax is not None) or (threshold is not None):
        # truncate the singular values
        chi = len(S)
        if chiMax is not None:
            chi = min(chi, chiMax)
        if threshold is not None:
            chi = min(chi, sum(np.cumsum(S**2)/sum(S**2) <= 1-threshold)+1)

        U = U[:, :chi]
        S = S[:chi]
        Vdg = Vdg[:chi, :]
    
    # normalize the singular values
    S = S / la.norm(S)

    return U, S, Vdg


```

````

Although we won't do it systematically in this course, MPS simulations should be repeated for different values of $\chi_\text{max}$ to determine where the simulations are converged. This is particularly important when studying ground states close to phase transitions, or when simulating unitary time evolution where the entanglement grows rapidly. In these cases, the required bond dimension can be very large, and it is important to know when the simulation is converged. In the context of TEBD, we will consider an example where we plot the entanglement entropy as a function of time for different values of $\chi_\text{max}$, which allows us to determine up to which time we can trust the simulation.
