# MPS Crash Course

In this course we will cover the basics of Matrix Product States (MPS) over four weeks. This will include: how MPS can be used to represent quantum many-body states, how to perform unitary time evolution using the TEBD algorithm, and how to find ground states using DMRG. 

This course is intended to be practical and will be based around writing your own code. You really won't get the most out of the course unless you try to write your own code following along with the lectures. However, I will provide a [GitHub repository with model code](https://github.com/MPS-crash-course/MPS-model-code). There will be a branch with model code for each week. The `main` branch (and branch called `week0-template`) provides a folder structure to start with, as well as yaml files for setting up a conda environment.

Although you will be writing your own code, I must stress that the goal is not to produce high-performance MPS code. Furthermore, the model code will provide the minimal code to cover the contents of this course. If you intend to use MPS (or tensor networks more generally) for your research, I would strongly recommend using a well-developed package. I can recommend [TeNPy](https://tenpy.readthedocs.io/en/latest/) for python, and [iTensor](https://itensor.org) for Julia, but there are many more. The goal of this course is to allow to use these packages without feeling like everything is a black box. At the very minimum, I hope that it demystifies those talks where the speaker very helpfully assumes that everyone is the audience is experts on MPS.

```{note}
This is the first time I am running this course, and I will be creating the content as I go along. I would be very grateful for any feedback on any aspect of the course so that I can improve it for future iterations.

I hope you find it useful!
```



```{tableofcontents}
```
