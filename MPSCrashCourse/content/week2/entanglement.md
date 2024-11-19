# Entanglement and Truncation

One of the central concepts in MPS is entanglement. Entanglement is a measure of how correlated two parts of a quantum system are. The entanglement spectrum (Schmidt values) also allows us to controllably truncate the MPS representation.

## The Schmidt decomposition

```{figure} images/schmidt.jpeg
---
name: fig:schmidt
width: 50%
align: center
---

???
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

???
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
        self.move_centre(site)

        ## YOUR CODE HERE ##

        return -np.sum(S**2 * np.log(S**2))


```

````


## Truncation

The Schmidt values provide us with a controlled way to truncate the MPS representation. The idea is to keep only the largest Schmidt values and discard the rest. This is equivalent to keeping only the most important correlations in the state.



````{admonition} Code: Entanglement entropy

???

```python
## file: src/svd.py

import numpy as np
import numpy.linalg as la

def svd_truncated(M, chiMax, threshold):

    U, S, Vdg = la.svd(M, full_matrices=False)
    
    if (chiMax is not None) and (threshold is not None):
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