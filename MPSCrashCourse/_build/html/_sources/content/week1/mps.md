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

```{figure} ../../images/simple_tensor_diagrams.jpeg
---
name: fig:simple_diagrams
width: 50%
align: center
---
Some simple objects as tensor network diagrams. (Left) A vector $v_i$. (Middle) A matrix $M_{ij}$. (Right) A rank-3 tensor $A_{ijk}$.
```

We can also represent a type of product between these tensors, which we call *contraction*. This is done diagramatically by connecting the legs of the two tensors. This corresponds to taking the product of the tensor elements and summing over the repeated index. For example, in {numref}`fig:contraction` we show the contraction of a matrix acting on a vector, and the contraction of two matrices, both in terms of matrix elements and as tensor network diagrams.


```{figure} ../../images/matrix_product_diagrams.jpeg
---
name: fig:contraction
width: 80%
align: center
---
(Left) A matrix acting on a vector written as a sum over indices and drawn as a tensor network diagram. (Right) Similar for the product of two matrices.
```

These kind of tensor contractions are the core part of the code we will write to manipulate MPS. 

We are then able to draw our MPS from Eq.{eq}`eq:mps` as a tensor network diagram, as shown in {numref}`fig:mps_diagram`. It consists of a line of rank-3 tensors. The legs that are not contracted are referred to as physical legs, and those that are contracted are virtual legs. The physical legs correspond to the indices of the tensor that are associated with the physical degrees of freedom of the system, in this case the spin states. The virtual legs correspond to the indices that are summed over in the contraction. 

```{figure} ../../images/mps_diagram.jpeg
---
name: fig:mps_diagram
width: 50%
align: center
---
A matrix product state on 5 sites as a tensor network diagram. The physical legs are the dangling vertical lines that are labelled, and the virtual legs are the horizontal lines. The end tensors are shown with a dashed line to indicate that there free indices are dimension 1. This allows us to write each tensor as a rank-3 numpy array, simplifying our code. Strictly speaking, the results of the contraction of this diagram is a $1\times 1$ matrix (which is equivalent to a scalar).
```



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

assert np.allclose(result, M @ v)  # check that the result matches

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

assert np.allclose(result1, A @ B)  # check that the result matches
assert np.allclose(result2, A @ B)  # check that the result matches

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

```{figure} ../../images/theta_example.jpeg
---
name: fig:theta_example
width: 50%
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

````

## MPS from a state vector

Now we have the diagramatic notation and the corresponding python functions at our disposal, let us return to the MPS representation. In particular, we want to explicitly relate the matrices in the MPS to the vector of probability amplitudes. The MPS can be drawn as the tensor network diagram shown in ???. Note the dashed virtual indices that are dangling on the end tensors. These are 1 dimensional indices and are simply included to make our code simpler (we only have to treat one shape of array). The result of this diagram when we specify all the physical indices is then not technically a number, but a $1\times 1$ matrix. 

We will now show how to convert a state vector into this MPS form. This will involve the use of numpy `reshape`, as well as singular value decomposition (SVD).

Let us start by thinking what we mean by treating $\psi_{i_1, i_2, \ldots, i_N}$ as a vector. One way to think of this is that $(i_1, i_2, \ldots, i_N)$ is a binary representation of the index of the vector. That is the probability amplitude corresponding to $|000\rangle$ is in the 0th position, $|001\rangle$ in the 1st position, and so on. Instead, we will no choose to separate out the first spin, and instead write the state as the matrix $\psi_{(i_1),(i_2 \ldots i_N)}$. That is, $|000\rangle$ corresponds to $\psi_{0, 0}$, $|001\rangle$ to $\psi_{0, 1}$, and $|100\rangle$ to $\psi_{1,0}$, and so on. This is a simple reshaping of the vector into a matrix, and can be done using `reshape` as follows in python:

```python

psi  # the vector of probability amplitudes

# Reshape the vector into a matrix
psi = psi.reshape(2, -1)

```

We specify that the matrix should have 2 rows, and then let python figure out how many columns are needed. Diagramatically this is shown in ???.



