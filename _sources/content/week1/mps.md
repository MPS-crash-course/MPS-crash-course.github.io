# Matrix Product States

Matrix product states are an efficient way of expressing and/or approximating many-body quantum states. They are particularly useful for studying one-dimensional quantum systems, but have applications in higher dimensions, and also for studying classical systems. In this course we will be focussing on the one-dimensional Heisenberg model, which is a model of interacting quantum spins. 

## Quantum states of many spins

A spin-1/2 degree of freedom is simply a two level quantum system. These two levels may correspond to the two possible spin states of a fundamental particle, such as an electron, or they could be effective degrees of freedom of a more complex system. 

```{note}
From now on we will use the terms "spin" and "spin-1/2 particle" interchangeably.
```

We typically label the two states as $|\uparrow\rangle$ and $|\downarrow\rangle$, or equivalently as $|0\rangle$ and $|1\rangle$. The most general state of a single spin can be written as a superposition of these two states, that is,

$$
|\psi\rangle = \psi_0 |0\rangle + \psi_1 |1\rangle,
$$

where $\psi_0$ and $\psi_1$ are complex numbers that satisfy $|\psi_0|^2 + |\psi_1|^2 = 1$. This is because the state $|\psi\rangle$ must be normalized, that is, $\langle \psi | \psi \rangle = 1$. These complex number $\psi_i$ are known as the *probability amplitudes* for the state $|\psi\rangle$. If we had two spins, then the most general state would be

$$
|\psi\rangle = \psi_{00} |00\rangle + \psi_{01} |01\rangle + \psi_{10} |10\rangle + \psi_{11} |11\rangle,
$$

where $\psi_{ij}$ satisfy the normalization condition $\sum_{i,j} |\psi_{ij}|^2 = 1$. Here we have introduced the shorthand notation $|ij\rangle = |i\rangle \otimes |j\rangle$ for the tensor product of the single spin states. Continuing in this fashion, the state of $N$ particles is given by

$$
|\psi\rangle = \sum_{i_1, i_2, \ldots, i_N} \psi_{i_1 i_2 \ldots i_N} |i_1 i_2 \ldots i_N\rangle,
$$ (eq:generalState)

where the sum is over all possible combinations of $i_1, i_2, \ldots, i_N \in \{\uparrow, \downarrow\}$. The number of probability amplitudes required to specify the state of $N$ particles grows exponentially with $N$, which is known as the curse of dimensionality. In exact numerical calculations, we typically store the complete list of probability amplitudes as a vector. For $N$ spins particles, this vector has $2^N$ complex numbers, which is becomes very quickly infeasible for large $N$.

## Matrix Product States

Matrix product states (MPS) provide a way of representing the state of a many-body quantum system in a more efficient way. The key idea is to write the state as a product of matrices, where each matrix corresponds to a single spin. The state of a single spin is then given by the product of these matrices. More explicitly, the probability amplitudes in {eq}`eq:generalState` can be rewritten in terms of a product of matrices $M^{[n] i_n}$ as

$$
\psi_{i_1 i_2 \ldots i_N} = M^{[1] i_1} M^{[2] i_2} \cdots M^{[N] i_N},
$$

where each $M^{[n] i_n}$ is a $\chi_n \times \chi_{n+1}$ matrix. We include the $[n]$ superscript to indicate that the matrix matrix is associated with the $n^{\text{th}}$ spin. The superscript $i_n$ indicates the state of the $n^{\text{th}}$ spin. That the probability amplitude (and hence the state) can be written in this way is perhaps not immediately obvious, but we will soon discuss how to find the matrices $M^{[n] i_n}$, and convert from a state vector to a MPS.

Instead of dealing with an MPS as a product of matrices, it is instead more useful to expose the matrix indices explicitly. That is $[M^{[n] i_n}]_{\alpha_n, \alpha_{n+1}} = M^{[n] i_n}_{\alpha_n, \alpha_{n+1}}$ where $\alpha_n$ and $\alpha_{n+1}$ are the indices of the matrix. The state of the system can then be written as

$$
\psi_{i_1 i_2 \ldots i_N} = \sum_{\alpha_1, \alpha_2, \ldots, \alpha_N} M^{[1]i_1}_{\alpha_0,\alpha_1} M^{[2]i_2}_{\alpha_2, \alpha_3}M^{[3]i_3}_{\alpha_3, \alpha_4} \cdots M^{[N]i_N}_{\alpha_{N},\alpha_{N+1}}.
$$ (eq:mps)

It is then more useful to treat the matrices $M^{[n]i_n}_{\alpha_n,\alpha_{n+1}}$, as rank-3 tensors, with indices $i_n, \alpha_n, \alpha_n$. Since these expressions can become cumbersome, we will also introduce a graphical notation to represent these tensors. This graphical notation is known as the tensor network diagram, and is a powerful tool for understanding and manipulating tensor networks. Let us introduce this notation before returning to this expression for the matrix product state.

````{admonition} Code: MPS Class

Since we know what the MPS representation of a state is, i.e., a collection of rank-3 tensors, one for each site, we can write a class that represents an MPS. We will start this class as follows:

```python
## file: src/mps.py

class MPS:
    """
    Matrix Product State class for 1D quantum systems of spin-1/2 particles.

    Attributes
    ----------
    L : Int 
        number of sites
    tensors : list of np.Array[ndim=3]
        list of tensors. Indices are (left, physical, right)
    """

    def __init__(self, L, tensors):
        self.L = L  # number of sites
        self.tensors = tensors  # list of tensors. Indices are (left, physical, right)
        

    def copy(self):
        return MPS(self.L, [tensor.copy() for tensor in self.tensors])
```

The class consists of two attributes: `L` which is the number of sites in the system, and `tensors` which is a list of rank-3 tensors. These tensors are numpy arrays. The `copy` method is a useful method to have, as we will often want to make a copy of an MPS. 

````


## Tensor network diagrams

Rather than constantly writing out tensors, their many indices, and the sums over these indices, we can introduce a simple diagramatic approach to represent these expressions, as shown in {numref}`fig:simple_diagrams`. In this notation, the tensor is represented by a shape (in this case a circle), and the outcoming lines (referred to as legs) represent the indices of the tensor. The number of legs is equal to the rank of the tensor. For example, the vector $v_i$ is represented by a circle with a single leg, the matrix $M_{ij}$ by two legs, and the rank-3 tensor $A_{ijk}$ by three legs. 

```{figure} images/simple_tensor_diagrams.jpeg
---
name: fig:simple_diagrams
width: 50%
align: center
---
Some simple objects as tensor network diagrams. (Left) A vector $v_i$. (Middle) A matrix $M_{ij}$. (Right) A rank-3 tensor $A_{ijk}$.
```

We can also represent a type of product between these tensors, which we call *contraction*. This is done diagramatically by connecting the legs of the two tensors. This corresponds to taking the product of the tensor elements and summing over the repeated index. For example, in {numref}`fig:contraction` we show the contraction of a matrix acting on a vector, and the contraction of two matrices, both in terms of matrix elements and as tensor network diagrams.


```{figure} images/matrix_product_diagrams.jpeg
---
name: fig:contraction
width: 80%
align: center
---
(Left) A matrix acting on a vector written as a sum over indices and drawn as a tensor network diagram. (Right) Similar for the product of two matrices.
```

These kind of tensor contractions are the core part of the code we will write to manipulate MPS. 

We are then able to draw our MPS from Eq.{eq}`eq:mps` as a tensor network diagram, as shown in {numref}`fig:mps_diagram`. It consists of a line of rank-3 tensors. The *physical legs* are those that are not connected, and correspond to the indices of the tensor that are associated with the physical degrees of freedom of the system, in this case the spin states. The *virtual legs* are the indices that are summed over in the contraction. 

```{figure} images/mps_diagram.jpeg
---
name: fig:mps_diagram
width: 50%
align: center
---
A matrix product state on 5 sites as a tensor network diagram. The physical legs are the dangling vertical lines that are labelled, and the virtual legs are the horizontal lines. The end tensors are shown with a dashed line to indicate that their free indices are dimension 1. This allows us to write each tensor as a rank-3 numpy array, simplifying our code. Strictly speaking, the results of the contraction of this diagram is a $1\times 1$ matrix (which is equivalent to a scalar).
```


(tensordot)=
## Numpy tensordot and transpose

In our MPS class, the tensors are given by rank-3 numpy arrays. In order to perform the contraction and keep track of the indices, we will use the numpy `tensordot` and `transpose` functions. The `tensordot` function is used to contract two tensors along specified axes. The `transpose` function is used to permute the axes of a tensor.

It is easiest to work through examples to understand how these functions work. Let use consider the examples in {numref}`fig:contraction`. Starting with the matrix acting on a vector, we can write the contraction as

```python
import numpy as np

# Define the tensors
v = np.array([1, 2])  # vector
M = np.array([[1, 2], [3, 4]])  # matrix

# Contract the matrix with the vector
result = np.tensordot(M, v, axes=([1], [0]))

# check that the result matches
assert np.allclose(result, M @ v), "matrix-vector multiplication failed!" 

```

The `tensordot` function takes three arguments. The first two are the arrays that we are contracting. The third argument specifies the axes that we are contracting over (specified as a tuple of lists). In this case we are contracting over the 1st of the matrix and 0th of the vector.

For the product of two matrices, let us do the contraction in two different ways by putting the matrices into `tensordot` in different orders.

```python

# Define the tensors
A = np.array([[1, 2], [3, 4]])  # matrix
B = np.array([[5, 6], [7, 8]])  # matrix

# Contract the matrices
result1 = np.tensordot(A, B, axes=([1], [0]))

# Contract the matrices in the opposite order
result2 = np.tensordot(B, A, axes=([0], [1]))
result2 = np.transpose(result2, (1, 0))  # transpose the result to match the order of the indices

# check that the result matches
assert np.allclose(result1, A @ B), "matrix-matrix method 1 failed!"  
assert np.allclose(result2, A @ B), "matrix-matrix method 2 failed!" 

```

When contracting using `tensordot`, order of indices will be a list of those remaining from the first array, followed by the remaining indices from the second. The two ways of contracting are equivalent to the following equations

$$
C_{ij} = \sum_k A_{ik} B_{kj}, \quad \text{and} \quad
\widetilde{C}_{ji} = \sum_k B_{kj} A_{ik}.
$$

In the second case we need to transpose the result to match the order of the indices. The `transpose` function allows us to perform a generalised transpose of tensors. The arguments are the array to transpose, followed by the new order of the axes.

Let us consider a final example using rank-3 tensors, which will be very similar to the type of contractions we will be performing in our MPS code.

```python

# Define the tensors
A = np.random.randn(4, 2, 4)  # (left, physical, right)
B = np.random.randn(4, 2, 4)  # (left, physical, right)

# Contract the tensors along the physical indices
theta = np.tensordot(A, B, axes=([1], [1]))  # (l1,p,r1) * (l2,p,r2) -> (l1,r1,l2,r2)
theta = np.transpose(theta, (0, 2, 1, 3))  # (l1,r1,l2,r2) -> (l1,l2,r1,r2)

```

Here we are computing the object

$$
\Theta_{ijkl} = \sum_p A_{i p k} B_{j p l},
$$

which we also show diagramatically in {numref}`fig:theta_example`. For more complex contractions it is good practice to add comments to your code keeping track of the indices. Here at the end we have transposed the result so that those on the left in {numref}`fig:theta_example` come before those on the right.

```{figure} images/theta_example.jpeg
---
name: fig:theta_example
width: 60%
align: center
---
The contraction of two rank-3 tensors $A$ and $B$ along the physical indices. The result is a rank-4 tensor $\Theta$. I have included coloured numbers showing the order of the indices as stored in the numpy arrays. To arrive at the final ordering for $\Theta$ we need to transpose the result.
```


````{admonition} Exercise: Contracting tensors

The `tensordot` and `transpose` functions can be quite confusing if you haven't dealt with tensors before, but are the main functions we will need throughout this course. I would therefore encourage you to play around with these to get a better feeling of how they work. The easiest way to do this is to work with matrices, where you can easily check the reaults by hand, or by printing the arrays.

Try the following exercises (all using square matrices complex matrices $A$ and $B$):

1) Compute $\text{Tr}(A B) = \sum_{ij} A_{ij} B_{ji}$ with $A$ in the first position and $B$ in the second.

2) Same as 1. but with the matrices appearing in the opposite order in `tensordot`.

3) Compute $A^\dagger B$ where $A^\dagger$ is the complex conjugate transpose of $A$. Compute this with $A^\dagger$ in the first position and $B$ in the second.

4) Compute $A^\dagger B$, but with $A^*$ (complex conjugate) in the second position and $B$ in the first.

I would recommend adding a file `tensordot.py` to your `exercises` directory and working through these exercises there. Make sure to add the line `from fix_pathing import root_dir`. You can then import your MPS class using `from src.mps import MPS`. 


````

## MPS from a state vector

Now we have the diagramatic notation and the corresponding python functions at our disposal, let us return to the MPS representation. In particular, we want to explicitly relate the matrices in the MPS to the vector of probability amplitudes. The MPS can be drawn as the tensor network diagram shown in {numref}`fig:mps_diagram`. Note the dashed virtual indices that are dangling on the end tensors. These are 1 dimensional indices and are simply included to make our code simpler (we only have to treat one shape of array). The result of this diagram when we specify all the physical indices is then not technically a number, but a $1\times 1$ matrix. 

We will now show how to convert a state vector into this MPS form. This will involve the use of numpy `reshape`, as well as singular value decomposition (SVD).

Let us start by thinking what we mean by treating $\psi_{i_1, i_2, \ldots, i_N}$ as a vector. One way to think of this is that $(i_1, i_2, \ldots, i_N)$ is a binary representation of the index of the vector. That is the probability amplitude corresponding to $|000\rangle$ is in the 0th position, $|001\rangle$ in the 1st position, and so on. Instead, we will now choose to separate out the first spin, and instead write the state as the matrix $\psi_{(i_1),(i_2 \ldots i_N)}$. That is, $|000\rangle$ corresponds to $\psi_{0, 0}$, $|001\rangle$ to $\psi_{0, 1}$, and $|100\rangle$ to $\psi_{1,0}$, and so on. This is a simple reshaping of the vector into a matrix, and can be done using `reshape` as follows in python:

```python

psi  # the vector of probability amplitudes

# Reshape the vector into a matrix
psi = psi.reshape(2, -1)

```

We specify that the matrix should have 2 rows, and then let python figure out how many columns are needed. Diagramatically this is shown in the first step of {numref}`fig:split_first_site`. We can then decompose this matrix using Singular Value Decomposition (SVD). The SVD of a matrix $A$ is given by

$$
A = U S V^\dagger,
$$

where $U$ and $V$ are unitary matrices, and $S$ is a diagonal matrix. The diagonal elements of $S$, typically labelled $\lambda_i$, are the singular values, which are non-negative, and are ordered from largest to smallest. In this context, the singular value decomposition is the same as the Schmidt decomposition, and the singular values are the Schmidt coefficients. Since the state we start from is pure and normalized, we have that $\sum_i \lambda_i^2 = 1$. These Schmidt values have an important physical meaning, which we will return to when we discuss truncation of the MPS next week. The SVD is shown diagramatically in {numref}`fig:split_first_site`, and can be computed in python as using `numpy.linalg` (or `scipy.linalg`) by

```python
import numpy.linalg as la

U, S, Vdg = la.svd(theta, full_matrices=False),
```
where we specify `full_matrices=False` since we are working with a rectangular matrix. If we have a $m \times n$ matrix, then this reduced SVD returns $U$ as $m \times k$, $S$ a length $k$ vector, and $V^\dagger$ as $k \times n$ matrix, where $k = \min(m,n)$.


```{figure} images/split_first_site.jpeg
---
name: fig:split_first_site
width: 60%
align: center
---

The first step of using Singular Value Decomposition to split a vector into an MPS. The vector is first reshaped into a matrix, and then decomposed into $U$, $S$, and $V^\dagger$. The matrix $U$ is reshapes into a rank-3 tensor and is the first tensor of the MPS. The matrices $S$ and $V^\dagger$ are combined to form the matrix $R$, which goes to the next step.
```

The matrix $U$ is then the first tensor in our MPS, i.e. $M^{[1]}$. For practical reasons, we want this tensor to be rank-3, so that we can deal with the tensors from each site in the same way. Therefore, we reshape this tensor again using `np.reshape(U, (1, 2, -1))`. The first index has dimension 1, and we represent this using a dashed line, as shown in {numref}`fig:split_first_site`. The second index has dimension 2 and this is the physical index, and the third index connects to the rest of the sites and is the virtual index. For the first site it will also be dimension 2. We then combine $S$ and $V^\dagger$ to form a tensor which we call $R$. 

```{figure} images/split_general_site.jpeg
---
name: fig:split_general_site
width: 70%
align: center
---

The general step for converting a state vector into an MPS. The tensor $R$ is reshaped and then decomposed using SVD. The resulting $U$ is reshaped into a rank-3 tensor, and the $S$ is multiplied into $V^\dagger$ to form the next $R$.
```


We then repeat this process down the chain, as shown in {numref}`fig:split_general_site`. We start by taking the tensor $R$, which should be a $\chi \times 2^m$ matrix, where $\chi$ is the number of Schmidt values from the previous step. We reshape this into a $2\chi \times 2^{m-1}$ tensor, and then perform the SVD. Effectively, this groups the left virtual index with the next physical site, and then the other dimension corresponds to the rest of the sites. We then perform the SVD, reshape $U$ to be $\chi \times 2 \times 2^{m-1}$, and multiply $S$ into $V^\dagger$ to form the next $R$. We repeat this process until we have tensors for all sites, as shown in {numref}`fig:state_to_mps`. When we get to the end of the chain, we reshape $R$ into a rank-3 tensor, and this is the final tensor in the MPS. 


```{figure} images/state_to_mps.jpeg
---
name: fig:state_to_mps
width: 100%
align: center
---

Tensor network diagrams for the process of converting a state vector to an MPS. The tensors of the MPS are split off from the vector by repeatedly reshaping and performing singular value decomposition. The final MPS is then a chain of tensors with bond dimension growing exponentially towards the centre of the chain.
```

This process of successively splitting the vector into tensors using SVD provides a proof-by-construction, that our MPS representation is equivalent to the original state vector. However, this process leads to a chain of tensor where the dimension of the virtual indices, known as the *bond dimension*, grows exponentially towards the centre of the chain. Hence, we have not actually gained anything. In practice, the power of MPS comes from truncating these tensors, restricting the bond dimension, and hence providing an approximation to the state vector. We will discuss this truncation process in more detail in the coming weeks.


````{admonition} Code: Extend the MPS Class

Since we have an algorithm to convert a state vector to an MPS, let us extend our MPS class to include this method. You should attempt this yourself, by following the structure below. A model solution can be found in the [GitHub repository](https://github.com/MPS-crash-course/MPS-model-code).

```python
## file: src/mps.py

class MPS:
    
    ## PREVIOUS CODE EXCLUDED ##

    @classmethod
    def fromVector(cls, vector):

        L = int(np.log2(vector.size))  # number of sites

        ## YOUR CODE HERE ##

        return cls(L, tensors)
```

We can add this as a `classmethod` to our MPS class. If you are not familiar with class methods, they are methods associated with the class type itself, not with an instance of the class. We would use this as follows: `psi = MPS.fromVector(vector)`. Here we are calling the method `fromVector`, on the abstract class type `MPS`, and it will return an instance of the class, which we call `psi`. Just like regular methods require `self` as the first argument, class methods require `cls`. 

````


## MPS to vector

We can also convert back from an MPS to a state vector. This is a simple process of contracting the tensors in the MPS, and then reshaping the result into a vector. We can do this by contracting the tensors from left to right, as shown in {numref}`fig:mps_to_state`. We start by contracting the first two tensors, then contract the result with the next tensor. Each time, we reshape the resulting rank-4 tensor into a rank-3 tensor by grouping the physical legs. We continue this process until we have contracted all the tensors, and then take the $(0,0)$ element of the resulting $1\times 1$ matrix.

```{figure} images/mps_to_state.jpeg
---
name: fig:mps_to_state
width: 100%
align: center
---

Tensor network diagrams for the process of converting an MPS to a state vector. The tensors of the MPS are contracted from left to right, and the result reshaped into a vector. 
```

````{admonition} Code: Extend the MPS Class

Let us also add the method to convert an MPS to a state vector.

```python
## file: src/mps.py

class MPS:
    
    ## PREVIOUS CODE EXCLUDED ##

    def toVector(self):

        ## YOUR CODE HERE ##

        return vector
```

We add this as a regular method to the MPS class. We can then call this method on an instance of the class, e.g. `vector = psi.toVector()`. 

````


## Product states

As we mentioned above, creating an MPS from a state vector is not very useful (except for testing our code, as we will below). The power of MPS comes from representing quantum states with a low, finite bond dimension. A simple example of states that can be exactly represented with minimal bond dimension are product states. We will use these both as a starting point for time evolution to simulate a quantum quench (week 2), but also as a starting point for the DMRG algorithm to find the ground state of the Heisenberg model (week 3).

A product state is a state which factorizes into a product of single site states. That is, the state of $N$ spins is given by

$$
|\psi\rangle = |\psi_1\rangle \otimes |\psi_2\rangle \otimes \cdots \otimes |\psi_N\rangle,
$$

where $|\psi_i\rangle$ is the state of the $i^{\text{th}}$ spin. The probability amplitudes for the many-body state $\psi_{i_1 i_2 \ldots i_N}$, can then be written as a product of the probability amplitudes for the single site states, i.e.,

$$
\psi_{i_1 i_2 \ldots i_N} = \psi_{i_1}^{[1]} \psi_{i_2}^{[2]} \cdots \psi_{i_N}^{[N]}.
$$

This product is actually already in the form of an MPS, with bond dimension 1. We simply need to take the two-dimensional vectors $\psi^{[i]}$ and reshape them into rank-3 tensors with dimensions $(1, 2, 1)$. We could easily add this to our MPS class (good exercise), but for now we will add a simplified version that creates a basis state, i.e., where each site is in a definite state, $\uparrow$ or $\downarrow$ ($0$ or $1$). 

````{admonition} Code: Product State

Let's add a method to the MPS class that creates a product state. We can specify the state as a list of $0$'s and $1$'s, where $0$ corresponds to $\uparrow$ and $1$ to $\downarrow$.

```python
## file: src/mps.py

class MPS:
    
    ## PREVIOUS CODE EXCLUDED ##

    @classmethod
    def productState(self, L, ):
        """
        Create a product state MPS on L sites in a given basis state on each site.

        Parameters
        ----------
        L : int
            Number of sites.
        state : list of ints
            List of basis states for each site. E.g. [0, 1, 0, 1] for a 4-site system.
        """

        tensors = []
        for s in state:
            assert s in [0, 1], "Basis states must be 0 or 1."

            if s == 0:
                tensors.append(np.array([1,0]).reshape((1,2,1)))
            else:
                tensors.append(np.array([0,1]).reshape((1,2,1)))

        return cls(L, tensors, 0)
```

````


## Summary

This week we have introduced the problem we want to solve and presented Matrix Product State methods as a way to solve it. We have introduced the basic structure of the MPS, and shown how to convert a state vector into an MPS and back again. In the process we introduced tensor network diagrams and the numpy functions `tensordot`, `transpose`, and `reshape`, which will be central to all the code we write in this course. It is important to understand these functions well, so I recommend playing about with these, and attempt the simple exercises in the Section {ref}`tensordot`.


````{admonition} Tests: MPS basics

It is also important that you test the MPS class and it's methods that you have written. Write the following tests.

1) Create a random complex vector of length $2^N$ (for $N<10$). Normalise this vector then convert it to an MPS and back to a vector. Check that the original vector and the final vector are the same.

2) Use the ED code (import using `from src.ed import *`) to compute the ground state of the Heisenberg model for $N<10$. Convert this state to an MPS and back to a vector. Check that the original vector and the final vector are the same.

3) Create the product state $|0000\rangle$ and convert this to a vector. Check that the vector is the same as the basis state $|0000\rangle$, it should have a 1 in the first element.

4) Create a basis state of your choice and convert to a vector. Check that the vector is the same as the basis state. The vector should have all zero elements except for a 1 in the position corresponding to the basis state. The binary representation of the basis state is the index of the 1 in the vector. e.g. $|0101\rangle$ should have a 1 in the 5th element of the vector.

I would recommend adding a file `mps_basics.py` to your `test` directory and working through these exercises there. Using $N<10$ should make this quick on most laptops. Make sure to add the line `from fix_pathing import root_dir`. You can then import your MPS class using `from src.mps import MPS`. 

````

Next week we will start by introducing a specific form of the MPS tensors called canonical form, which will allow us to perform more efficient calculations, and is important for correctly truncating our tensors. We will then implement the Time Evolving Block Decimation (TEBD) algorithm, which is a way to simulate the time evolution of a quantum system using MPS. 

