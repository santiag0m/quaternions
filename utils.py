from math import sin, cos, acos, radians
import quaternions.Quaternion as qt


def gyro2quaternion(g_xyz, fs=60):
    """
    Transforms gyroscope measurements (angular velocity: deg/s) to a rotation quaternion
    :param g_xyz: Gyro measurements (degrees)
    :param fs: Sensor sampling frequency
    :return: Quaternion object representing the rotation
    """
    assert len(g_xyz) == 3
    dt = fs ** -1
    g_xyz = [radians(w) for w in g_xyz]
    magnitude = qt.vector_norm(g_xyz)
    a = dt * magnitude * 0.5
    if magnitude > 0:
        x, y, z = qt.vector_scale(g_xyz, magnitude ** -1)
    else:
        x, y, z = 0
    s = sin(a)
    w = cos(a)
    x *= s
    y *= s
    z *= s
    return qt.Quaternion(w, x, y, z)


def tiltaa(q_a, q_t):
    """
    Returns the tilt correction quaternion for accelerometer data in angle-axis format
    :param q_a: Accelerometer vector quaternion
    :param q_t: Rotation quaternion from body to inertial frame
    :return: Tilt correction axis and angle
    """
    q_a = q_a.rotate(q_t)
    q_a.normalize()
    w, x, y, z = q_a.components()
    # phi = acos(y)  # Ideal component (Global reference)
    phi = acos(x)  # Ideal component (Initial reference)
    # v_tilt = cross_product((x, y, z), (0, 1, 0))  # Gravity vector (Global reference)
    v_tilt = qt.cross_product((x, y, z), (1, 0, 0))  # Gravity vector (Initial reference)
    return phi, v_tilt


def omegaconstruction(w):
    """
    Method to construct a symplectic-like matrix from a 3D vector. Use case for representing
    quaternion multiplication for update from gyro data (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5621018/)

    Om(w) = | 0  -wT |
            | w  [wX]|

    :param w: 3D vector to construct the matrix or a quaternion
    :return: Resulting block matrix
    """
    if isinstance(w, qt.Quaternion):
        w, x, y, z = w.components()
        v = [x, y, z]
        out = [[w] + list(qt.vector_scale(v, -1))]
        ss = qt.identity(3) * w
        ss += qt.Matrix(skewsymmetric(v))
        ss = ss.values
        ss_expanded = [[v[i]] + ess for i, ess in enumerate(ss)]
        out += ss_expanded
        return out
    else:
        assert len(w) >= 3, "The omega matrix block construction is only implemented for 3D vectors and quaternions"
        assert len(w) <= 4, "The omega matrix block construction is only implemented for 3D vectors and quaternions"
        if len(w) == 3:
            omegaconstruction(qt.Quaternion(0, *w))
        else:
            omegaconstruction(qt.Quaternion(*w))


def skewsymmetric(w):
    """
    Method to construct a skew-symmetric matrix from a 3D vector

    [wX] =  |  0  -w_z  w_y |
            | w_z   0  -w_x |
            |-w_y  w_x   0  |

    :param w: 3D vector to construct the matrix
    :return: Resulting skew-symmetric matrix
    """
    assert len(w) == 3, "The skew symmetric matrix is only implemented for 3D vectors"
    out = [[0, -1 * w[2], w[1]],
           [w[2], 0, -1 * w[0]],
           [-1 * w[1], w[0], 0]]
    return out


