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
$$

It is then more useful to treat the matrices $M^{[n]i_n}_{\alpha_n,\alpha_{n+1}}$, as rank-3 tensors, with indices $i_n, \alpha_n, \alpha_n$. Since these expressions can become cumbersome, we will also introduce a graphical notation to represent these tensors. This graphical notation is known as the tensor network diagram, and is a powerful tool for understanding and manipulating tensor networks. Let us introduce this notation before returning to this expression for the matrix product state.

## Tensor network diagrams

