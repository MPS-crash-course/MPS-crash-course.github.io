# Canonical Form and Truncation

Last week, we showed that matrix product state allows us to represent quantum many-body states in terms of a product of matrices (or as a chain of tensors). This representation is not unique, and we are free to change the individual tensors in the MPS without changing the state it represents. This freedom is known as gauge freedom. For instance, on each virtual bond, we can apply the identity $1 = X^{-1} X$, where $X$ is an invertible matrix. By contracting $X^{-1}$ and $X$ with the tensors on the left and right of the bond, we can change the tensors without changing the state, as shown in {numref}`fig:gauge`. We have this freedom on each virtual bond of the MPS

```{figure} images/gauge.jpeg
---
name: fig:gauge
width: 90%
align: center
---

???
```

Given this freedom in the MPS representation, there are particular forms for the tensors that have special useful properties that simplify calculations. These are called *canonical forms*. The most commonly used canonical forms are the left and right canonical forms. We will label tensors $A^{[n]}$ if they are in left canonical form, and $B^{[n]}$ if they are in right canonical form. Left canonical form is defined by the conditions shown in tensor network notation in {numref}`fig:canonical`. That is, contracting left canonical form tensors from the left results in the identity. Additionally, the diagonal matrix of Schmidt values is the unique right eigenvector of the left-canonical transfer matrix. This second condition is not one that we will make use of in this course, but it is a useful property of the left canonical form. The right canonical form is defined similarly, but with the conditions reversed.

```{figure} images/canonical.jpeg
---
name: fig:canonical
width: 75%
align: center
---

???
```

A state can be written purely in left or right canonical form, or we can use a mixed canonical form, as shown in {numref}`fig:mixed`. On the bond between the left and right canonical tensors, we have the Schmidt values (singular values of the Schmidt decomposition). It is common to keep these Schmidt values on the bond between the left and right canonical tensors, as this allows us to easily compute the entanglement entropy of the state. However, to simplify our code, we will instead introduce a central tensor, which we will call $C^{[n]}$. This tensor is equivalent to a left canonical tensor with the Schmidt values contracted from the right, and also to a right canonical tensor with the Schmidt values contracted from the left. We do this so that our state is still simply represented by a list of tensors all of the same rank. In order to do this, we need to keep track of the location of this central tensor. We add an attribute to our class called `centre`, which is the index of the central tensor. If the central tensor is at index $n$, then the left canonical tensors are at indices $1$ to $n-1$, the central tensor is at index $n$, and the right canonical tensors are at indices $n+1$ to $N$. We will later add a method to our class that moves the central tensor to a specified index.

```{figure} images/mixed.jpeg
---
name: fig:mixed
width: 100%
align: center
---

???
```

## Computing expectation values

So, how does canonical form help us? Let us consider the computation of a local expectation value, e.g., $\langle \psi | \sigma^z_n | \psi \rangle$. Expectation value can be written as the tensor network diagram shown in {numref}`fig:expectation`. If we did not use canonical form, then we would need to perform the contraction site-by-site all the way along the length of the chain. However, if we use a mixed canonical form where the centre is on the site of the operator, the contraction to the left trivially gives us identity, and the same to the right. Therefore, we are left to contract only three tensors on the site of the operator, dramatically reducing the computation cost.


```{figure} images/expectation.jpeg
---
name: fig:expectation
width: 90%
align: center
---

???
```