# quaternions

This repository contains auxiliary functions for quaternion based motion tracking written entirely in base python.

The base object 'Quaternion' can be imported directly using:

```
from quaternions.Quaternion import Quaternion

my_quaternion = Quaternion(q_0, q_1, q_2, q_3)

```

However, some auxiliary functions and the 'Matrix' class won't be available. So it is recommended to import the library in the following way:

```
import quaternions.Quaternion as qtn

my_quaternion = qt.Quaternion(q_0, q_1, q_2, q_3)

my_matrix = qt.Matrix([[1, 2], [3, 4]])

```

Matrix and Quaternion usage is designed to be as similar as numpy array usage as possible.
