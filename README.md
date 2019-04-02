# quaternions

This repository contains auxiliary functions for quaternion based motion tracking written entirely in base python.

## Basics

The main object `Quaternion` can be imported directly using:

```python
from quaternions.Quaternion import Quaternion

my_quaternion = Quaternion(q_0, q_1, q_2, q_3)

```

However, some auxiliary functions and the `Matrix` class won't be available. So it is recommended to import the library in the following way:

```python
import quaternions.Quaternion as qtn

my_quaternion = qtn.Quaternion(q_0, q_1, q_2, q_3)

my_matrix = qtn.Matrix([[1, 2], [3, 4]])

```

## Rotations

### From and To Tait-Bryan Angles

A rotation quaternion can be created from the Tait-Bryan angles using:

```python
from math import radians

yaw = radians(90)
pitch = radians(45)
roll = radians(30)

rot_quaternion = qtn.tb2quaternion(yaw, pitch, roll)
```
And a vector quaternion can be rotated by:

```python
new_quaternion = my_quaternion.rotate(rot_quaternion)
```
Or transform the rotation quaternion back to Tait-Bryan angles:

```python
yaw, pitch, roll = rot_quaternion.totb
```

It must be noted that the coordinate axis used is the standard coordinate axis (x - forward, y - left, z - up) and rotation order is Z (yaw) -> Y (pitch) -> X (roll).

### From and To Angle-Axis notation

Rotation quaternions can also be created using the angle-axis notation:

```python
angle = radians(90)
axis = (1, 1, 1)
rot_quaternion = qtn.aa2quaternion(angle, axis)
```

And back:

```python
angle, axis = rot_quaternion.aa
```

The resulting quaternion will aply a rotation of `angle` radians around the 3D `axis` vector.

## Notes

Remember that all vector quaternions should have real component (`q_0`) equal to 0.

Remember that all rotation quaternions must be normalized. This can be done using the in-place operation `normalize`:

```python
rot_quaternion.normalize()
```

Matrix and Quaternion usage is designed to be as similar as numpy array usage as possible, however 2D indexing has not been implemented yet.

