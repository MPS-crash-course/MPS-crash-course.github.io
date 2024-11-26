# Matrix Product Operators

This week we will introduce the Density Matrix Renormalization Group (DMRG) algorithm to find the ground state of a 1D quantum many-body system. It may seem that finding ground states would be the first task to tackle before time evolution, especially for the problem we set out. However, the DMRG algorithm is conceptually more complex than TEBD. The main reason for this is that the modern formulation of DMRG is based on the Matrix Product Operator (MPO) representation of the Hamiltonian, and so in this section we will introduce MPOs.

```{figure} images/mpo.jpeg
---
name: fig:mpo
width: 65%
align: center
---

The MPO representation of an operator, such as the Hamiltonian. The operator is represented as a product of rank-4 tensors.
```

Similar to MPS, MPOs are a way to represent the elements of an operator as a product of matrices. It will again be more useful to consider the full operator at the product of rank-4 tensors. The MPO representation of an operator is shown in {numref}`fig:mpo`. As with MPS, we also include a dimension 1 index on the first and last tensors so that all tensors have the same rank. In our code, we will choose the labelling of the legs as shown in {numref}`fig:mpo_tensor`. 


```{figure} images/mpo_tensor.jpeg
---
name: fig:mpo_tensor
width: 33%
align: center
---

Labeling the indices of the rank-4 tensors in the MPO representation.
```

```{note}
Unlike MPS, we don't have a canonical form for MPOs. 
```


````{admonition} Code: MPO Class

Let us start writing the MPO class. It will consist of a set of tensors, and the number of sites. For now let us simply writing the initialization.

```python
## file: src/mpo.py

import numpy as np

class MPO:
    """
    Matrix Product State class for 1D quantum systems of spin-1/2 particles.
    """

    def __init__(self, L, tensors):
        self.L = L  # number of sites
        self.tensors = tensors  # list of tensors. Indices are (left, p_out, p_in, right)

```

````



## MPO representation of the Hamiltonian

```{figure} images/FSM.jpeg
---
name: fig:FSM
width: 60%
align: center
---

Finite State Machine representation of the Heisenberg Hamiltonian. The MPO tensors can be easily extracted from this diagram.
```

All local Hamiltonians can be written *exactly* as an MPO with a finite bond dimension. In the case of the Heisenberg model, we can write it as an MPO with bond dimension 5. To construct the MPO representation of the Hamiltonian, we use a Finite State Machine (FSM) representation of the Hamiltonian, as shown in {numref}`fig:FSM`. See Ref. {cite}`Crosswhite2008` for more details of this construction. Ultimately, we end up with the following MPO tensors for the Heisenberg model:

$$
W^{[n]} = \left(\begin{matrix}
1 & X & Y & Z & 0 \\
0 & 0 & 0 & 0 & X \\
0 & 0 & 0 & 0 & Y \\
0 & 0 & 0 & 0 & Z \\
0 & 0 & 0 & 0 & 1
\end{matrix}\right)
$$ (eq:mpo_heisenberg)

where $X$, $Y$, and $Z$ are the Pauli matrices. Here we have written the MPO tensor as a matrix of matrices. The outer matrix correspond to the virtual indices, and the inner matrices correspond to the physical indices. More explicitly, we have, e.g., $W^{[n]}_{0,i,j,0} = 1_{i,j}$ and $W^{[n]}_{0,i,j,1} = X_{i,j}$.


````{admonition} Code: Create Hamiltonian MPO

We can now add the class method to create the MPO for the Heisenberg model.

```python
## file: src/mpo.py

class MPO:
    
    ## PREVIOUS CODE OMITTED ##

    @classmethod
    def Hamiltonian(cls, L):
        """
        Construct the MPO for the 1D Heisenberg Hamiltonian of length L.
        """
        # Define the Pauli matrices
        identity = np.eye(2)
        X = np.array([[0, 1], [1, 0]])
        Y = np.array([[0, -1j], [1j, 0]])
        Z = np.array([[1, 0], [0, -1]])

        W = np.zeros((5,2,2,5), dtype=complex)
        W[0, :, :, 0] = identity
        W[0, :, :, 1] = X
        W[0, :, :, 2] = Y
        W[0, :, :, 3] = Z
        W[1, :, :, 4] = X
        W[2, :, :, 4] = Y
        W[3, :, :, 4] = Z
        W[4, :, :, 4] = identity

        # Construct the Heisenberg Hamiltonian
        tensors = [W.copy() for _ in range(L)]

        tensors[0] = tensors[0][0, :, :, :].reshape(1, 2, 2, 5)
        tensors[-1] = tensors[-1][:, :, :, 4].reshape(5, 2, 2, 1)

        return cls(L, tensors)


```

````



## Expectation values


```{figure} images/mpo_expectation.jpeg
---
name: fig:mpo_expectation
width: 80%
align: center
---

Computation of the expectation value of the Hamiltonian using the MPS and MPO representations. To contract the diagram we first contract the slices and convert them to matrices. The full contraction is then the matrix product of all the slices.
```



---

## References

```{bibliography}
:filter: docname in docnames
```