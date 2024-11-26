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
1 & S^x & S^y & S^z & 0 \\
0 & 0 & 0 & 0 & S^x \\
0 & 0 & 0 & 0 & S^y \\
0 & 0 & 0 & 0 & S^z \\
0 & 0 & 0 & 0 & 1
\end{matrix}\right)
$$ (eq:mpo_heisenberg)

where $S^\alpha = \frac{1}{2}\sigma^\alpha$ are the spin operators. Here we have written the MPO tensor as a matrix of matrices. The outer matrix correspond to the virtual indices, and the inner matrices correspond to the physical indices. More explicitly, we have, e.g., $W^{[n]}_{0,i,j,0} = 1_{i,j}$ and $W^{[n]}_{0,i,j,1} = X_{i,j}$.


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
        # Define the spin matrices
        identity = np.eye(2)
        Sx = 1/2*np.array([[0, 1], [1, 0]])
        Sy = 1/2*np.array([[0, -1j], [1j, 0]])
        Sz = 1/2*np.array([[1, 0], [0, -1]])

        W = np.zeros((5,2,2,5), dtype=complex)
        W[0, :, :, 0] = identity
        W[0, :, :, 1] = Sx
        W[0, :, :, 2] = Sy
        W[0, :, :, 3] = Sz
        W[1, :, :, 4] = Sx
        W[2, :, :, 4] = Sy
        W[3, :, :, 4] = Sz
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

We will also want to compute the expectation value of the Hamiltonian in MPO form. This is done by contracting the MPO with the MPS and its complex conjugate, as shown in {numref}`fig:mpo_expectation`. In this diagram, we can choose to perform the contraction in many different orders. We will perform the contraction by first contracting the "slices" of the MPO and MPS (i.e. the transfer matrix). When reshaped as matrices, the full contraction is then the matrix product of all the slices. We will also need to compute these slices for the DMRG algorithm. For these reasons, we will add a method to our class to construct the slice at a given site. We can then use that to compute the expectation value of the Hamiltonian.

````{admonition} Code: MPO Expectation Value

Let us add two new methods to the MPO class. The first computed the slice (transfer matrix) at site $i$ and the second computes the expectation value of the Hamiltonian. The ordering of the indices and contractions that I used are shown in Fig. {numref}`fig:slice_indexing`.

```{figure} images/slice_indexing.jpeg
---
name: fig:slice_indexing
width: 70%
align: center
---

Contraction of MPS and MPO tensors to construct the slice (transfer matrix) at site $i$. The red indices correspond to the first input to `tensordot`, and green are the second. The final step is a transpose of the indices.
```


```python
## file: src/mpo.py

class MPO:
    
    ## PREVIOUS CODE OMITTED ##

    def get_slice(self, psi, i):
        ## YOUR CODE HERE ##

        # return the slice as a rank-6 tensor
        return M


    def expectation(self, psi):
        """
        Compute the (real) expectation value of the (Hermitian) MPS with respect to a state.
        """
        assert psi.L == self.L, "State size does not match MPS size."

        ## YOUR CODE HERE ##

        return np.real(overlap[0, 0])   

```
````



````{admonition} Tests: MPO

Let us add a simple test for the MPO class. We will compute the expectation value of the Hamiltonian for a random state and compare it to the exact value computed using our ED code. Your test can use the following steps:

1) Create a random normalized state vector.

2) Convert the state vector to an MPS.

3) Construct the exact Hamiltonian matrix using the ED code. Use this matrix to compute the expectation value with respect to the state vector from step 1.

4) Use the MPO class method to construct the Hamiltonian MPO. Use the method on this instance to compute the expectation value with respect MPS from step 2.

5) Compare the two expectation values.

I would recommend adding a file `mpo.py` to your `test` directory and adding this test. Using $L<10$ should make this quick on most laptops. Make sure to add the line `from fix_pathing import root_dir`. 

````





---

## References

```{bibliography}
:filter: docname in docnames
```