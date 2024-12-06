# Beyond This Course

I hope you have enjoyed this crash course into the basics of MPS. Naturally, there is much more to learn about MPS, both in terms of numerical methods, but also as an analytical tool, and in the more general framework of tensor networks. Here are just a few resources to help you continue your journey:

- {cite}`Hauschild2018` is a great set of lecture notes that covers many of the topics we have discussed in this course, and more. In particular, it also discussed infinite generalisations of MPS, TEBD and DMRG. These lecture notes were the main inspiration for this course.

- {cite}`Schollwock2011,Orus2014,Cirac2021`. These reviews are a great starting point for learning more about tensor networks and matrix product states. They cover a wide range of topics and are written by experts in the field.

- [TeNPy](https://tenpy.readthedocs.io/en/latest/) is the documentation for the `tenpy` library, which is a Python library for tensor network calculations. 

- [iTensor](https://itensor.org/) is a Julia library for tensor network calculations. 


## Developing your code further

The purpose of this course was to walk you through the basics of MPS by writing your own code. I feel like this is the best way to learn the topic if you are planning to use them as a numerical method. Beyond this, I would definitely recommend the use of a library like `tenpy` or `iTensor` for more advanced calculations. However, if you are interested in developing your own code further, here are a few ideas:

- **Iterative Eigensolver**: The biggest bottleneck in the DMRG algorithm is solving the local eigenproblem. The code we wrote can be sped up significantly by using an iterative eigensolver such as the Davidson algorithm. This is a relatively simple algorithm. By using very few iterations in the eigensolver, you can afford to do more sweeps in the DMRG algorithm, to get a comparable accuracy with sufficiently smaller computational cost.

- **Add convergence criteria**: We manually controlled how many sweeps to do in the DMRG algorithm, but it is possible that we could have stopped earlier. You could add a convergence criterion to the code, such as checking when the change in the energy between sweeps is below a certain threshold.

- **Perform convergence tests**: We didn't test how the energy converged with the bond dimension. You could add this to your code to see how the energy changes as you increase the bond dimension for DMRG. For the time evolution you should track the entanglement entropy and plot this as a function of time for different bond dimensions.

- **Implement different Hamiltonian MPOs**: We constructed the MPO for the Antiferromagnetic Heisenberg chain. However, this construction can easily be generalised to other models, e.g. the XXZ model or the Ising model with field. This will open up the possibility of finding the ground states of a wider range of models.

- **Implement a more general TEBD step**: Our TEBD algorithm was hard coded for the AFH chain. You can generalise this for other Hamiltonians with two-body nearest neighbour Hamiltonians.

- **Implement other time evolution methods**: The TEBD algorithm is just one way of performing time evolution. As we introduced it, it is limited to local nearest neighbour Hamiltonians. You could implement other methods, such as those based on MPOs, which allow you to simulate more general Hamiltonians, including long-range interactions. See e.g. {cite}`Zaletel2015`.


## Feedback on the course

I really hope you have enjoyed this course and that it will be useful for your future studies. This is the first time I have offered this course, but I now see that it could be a valuable resources for others and so I may turn this into a MPAGS course in the future. I would love to hear your feedback on the course, and any suggestions you have for improvements. Please feel free to email me at [adam.gammon-smith@nottingham.ac.uk](mailto:adam.gammon-smith@nottingham.ac.uk) with any feedback you have.

---

## References

```{bibliography}
:filter: docname in docnames
```