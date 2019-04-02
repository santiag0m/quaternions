# quaternions

This repository contains auxiliary functions for quaternion based motion tracking written entirely in base python.

The main object 'Quaternion' can be imported directly using:

```python
from quaternions.Quaternion import Quaternion

my_quaternion = Quaternion(q_0, q_1, q_2, q_3)

```

However, some auxiliary functions and the 'Matrix' class won't be available. So it is recommended to import the library in the following way:

```python
import quaternions.Quaternion as qtn

my_quaternion = qtn.Quaternion(q_0, q_1, q_2, q_3)

my_matrix = qtn.Matrix([[1, 2], [3, 4]])

```

A rotation quaternion can be created from the Tait-Bryan angles using:

```python
from math import radians

yaw = radians(90)
pitch = radians(45)
roll = radians(30)

rot_quaternion = qtn.tb2quaternion(yaw, pitch, roll)
```
And a vector quaternion can be rotated by:

```
new_quaternion = my_quaternion.rotate(rot_quaternion)
```

It must be noted that the coordinate axis used is the standard coordinate axis (x - forward, y - left, z - up) and rotation order is Z (yaw) -> Y (pitch) -> X (roll).

Matrix and Quaternion usage is designed to be as similar as numpy array usage as possible.
