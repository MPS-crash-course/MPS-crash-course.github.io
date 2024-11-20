# Canonical Form

Last week, we showed that matrix product state allows us to represent quantum many-body states in terms of a product of matrices (or as a chain of tensors). This representation is not unique, and we are free to change the individual tensors in the MPS without changing the state it represents. This freedom is known as gauge freedom. For instance, on each virtual bond, we can apply the identity $1 = X^{-1} X$, where $X$ is an invertible matrix. By contracting $X^{-1}$ and $X$ with the tensors on the left and right of the bond, we can change the tensors without changing the state, as shown in {numref}`fig:gauge`. We have this freedom on each virtual bond of the MPS

```{figure} images/gauge.jpeg
---
name: fig:gauge
width: 90%
align: center
---

Gauge freedom in defining the tensors of an MPS. By introducing an invertible matrix (and its inverse) on the virtual bonds of the MPS, we can change the tensors without changing the state.
```

Given this freedom in the MPS representation, there are particular forms for the tensors that have special useful properties that simplify calculations. These are called *canonical forms*. The most commonly used canonical forms are the left and right canonical forms. We will label tensors $A^{[n]}$ if they are in left canonical form, and $B^{[n]}$ if they are in right canonical form. Left canonical form is defined by the conditions shown in tensor network notation in {numref}`fig:canonical`. That is, contracting left canonical form tensors from the left results in the identity. Additionally, the diagonal matrix of Schmidt values is the unique right eigenvector of the left-canonical transfer matrix. This second condition is not one that we will make use of in this course, but it is a useful property of the left canonical form. The right canonical form is defined similarly, but with the conditions reversed.

```{figure} images/canonical.jpeg
---
name: fig:canonical
width: 75%
align: center
---

Left and right canonical forms. Contracting the left virtual and physical indices of left canonical tensors with their conjugate results in the identity. The square of the Schmidt values are the right eigenvector of the left canonical form transfer matrix. Similarly for right canonical, but with in other direction.
```

A state can be written purely in left or right canonical form, or we can use a mixed canonical form, as shown in {numref}`fig:mixed`. On the bond between the left and right canonical tensors, we have the Schmidt values (singular values of the Schmidt decomposition). It is common to keep these Schmidt values on the bond between the left and right canonical tensors, as this allows us to easily compute the entanglement entropy of the state. However, to simplify our code, we will instead introduce a central tensor, which we will call $C^{[n]}$. This tensor is equivalent to a left canonical tensor with the Schmidt values contracted from the right, and also to a right canonical tensor with the Schmidt values contracted from the left. We do this so that our state is still simply represented by a list of tensors all of the same rank. In order to do this, we need to keep track of the location of this central tensor. We add an attribute to our class called `centre`, which is the index of the central tensor. If the central tensor is at index $n$, then the left canonical tensors are at indices $1$ to $n-1$, the central tensor is at index $n$, and the right canonical tensors are at indices $n+1$ to $N$. We will later add a method to our class that moves the central tensor to a specified index. 

```{figure} images/mixed.jpeg
---
name: fig:mixed
width: 100%
align: center
---

Representing an MPS in mixed canonical form. We choose to incorporate the Schmidt values into a central tensor, which we will call $C^{[n]}$. This tensor is equivalent to a left canonical tensor with the Schmidt values contracted from the right, and also to a right canonical tensor with the Schmidt values contracted from the left. We will keep track of the location of this central tensor with an attribute `centre` in our MPS class.
```

```{note}

Our code that converted a state vector into an MPS resulted in tensors in the left canonical form. This is because we started from the left and successively performed SVD, leaving behind left canonical tensors. Therefore, the centre of the MPS is at the last site of the chain. 

The product state that we constructed is special. All of its tensors are simultaneously in left and right canonical form. This is because there is a single Schmidt value equal to 1 on each bond. We can therefore arbitrarily choose which site we call the centre.

```

````{admonition} Code: Add centre to MPS class

In order to keep track of the centre site, we need to add an attribute to our MPS class. We will also need to update our previous methods, since the class initialization will now have an extra argument. We therefore have to make the following changes.

```python

class MPS:

    def __init__(self, L, tensors, centre):
        self.L = L  # number of sites
        self.tensors = tensors  # list of tensors. Indices are (left, physical, right)
        self.centre = centre  # position of the orthogonality centre
        
    def copy(self):
        return MPS(self.L, [tensor.copy() for tensor in self.tensors], self.centre)

    @classmethod
    def fromVector(cls, vector):
        ## YOUR CODE HERE

        return cls(L, tensors, L-1)

    @classmethod
    def productState(cls, L, state):
        ## YOUR CODE HERE

        return cls(L, tensors, 0)

```

This will be the final set of attributes for the MPS class.

````



## Computing expectation values

So, how does canonical form help us? Let us consider the computation of a local expectation value, e.g., $\langle \psi | \sigma^z_n | \psi \rangle$. Expectation value can be written as the tensor network diagram shown in {numref}`fig:expectation`. If we did not use canonical form, then we would need to perform the contraction site-by-site all the way along the length of the chain. However, if we use a mixed canonical form where the centre is on the site of the operator, the contraction to the left trivially gives us identity, and the same to the right. Therefore, we are left to contract only three tensors on the site of the operator, dramatically reducing the computation cost.


```{figure} images/expectation.jpeg
---
name: fig:expectation
width: 90%
align: center
---

The computation of the expectation value for a single site operator. Using canonical form, the contraction of the two copies of the MPS reduces to a simple tensor contraction on the site of the operator.
```



## Moving the centre

Since the computation of expectation values is most efficient when the operator acts on the same site as the position of the central tensor, we need to be able to move the centre to a desired site. This will also be necessary when computing the entanglement entropy. To do this we will use singular value decomposition to separate the Schmidt values from the central tensor, allowing us to contract them to the left or right. For example, in {numref}`fig:move_centre`, we show the process of moving the centre to the right. We start by contracting the central tensor with the tensor to its right (which should be in right canonical form). We then perform an SVD on the resulting tensor. By contracting the Schmidt values with the right tensor, we end up with a left canonical form tensor followed by the central tensor, and hence have moved the centre by one site to the right. An analogous process can be used to move the centre to the left.

```{figure} images/move_centre.jpeg
---
name: fig:move_centre
width: 100%
align: center
---

The computation of the expectation value for a single site operator. Using canonical form, the contraction of the two copies of the MPS reduces to a simple tensor contraction on the site of the operator.
```

Note that when we contract the tensors to move the centre, we will want to keep track of the bond dimension $\chi$ of the leg that we contract. This is because the SVD may result in tensors with larger bond dimension. However, we know that only the first $\chi$ Schmidt values will be non-zero, so we should truncate the Schmidt values and the tensors to keep the bond dimension fixed. This truncation does not incur any approximation. Without doing this, our bond dimension would grow out of control. We will cover truncation for approximation purposes in the next section.



````{admonition} Code: Move centre and compute expectation value

We can now add new methods to our MPS class that move the centre to the left or right, or to any site we choose. We can also add a method that computes the expectation value of a local operator acting on a site. 

```python

class MPS:

    ## Previous code omitted

    def move_centre_left(self):
        ## YOUR CODE HERE
        ## Remember to truncate the bond dimension to match the original bond dimension

        self.centre -= 1

    def move_centre_right(self):
        ## YOUR CODE HERE
        ## Remember to truncate the bond dimension to match the original bond dimension

        self.centre += 1

    def move_centre_to(self, i):
        """
        Move the orthogonality centre to site i.
        """

        while self.centre > i:
            self.move_centre_left()

        while self.centre < i:
            self.move_centre_right()

    def expectation(self, O, site):
        """
        Compute the expectation value of an operator at a given site.

        Parameters
        ----------
        O : np.ndarray (shape=(2,2))
            Operator acting on a single site.
        site : int
            Site at which to compute the expectation value.
        """

        self.move_centre_to(site)

        ## YOUR CODE HERE

        return expectation 

```

When you move the centre you want to catch the case where the centre is already at the end of the chain. Otherwise the while loop will continue indefinitely. You may also replace the while loops for safety.

Remember to truncate the bond dimension of the tensor after the SVD to match the original bond dimension. This increase in bond dimension does not contain additional information, and so we can safely truncate it without approximation.

````



