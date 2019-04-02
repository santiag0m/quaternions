from math import atan2, asin, acos, sin, cos


class Matrix:
    def __init__(self, array):
        """
        A 2D Matrix object
        :param array: List-of-lists or tuple-of-tuples values
        """
        assert isinstance(array, (list, tuple)), "Input should be a list-of-lists (or tuples)"
        assert isinstance(array[0], (list, tuple)), "Input should be a list-of-lists (or tuples)"
        c = len(array[0])
        for r in array:
            assert len(r) == c, "All columns should have the same number of elements"
        self.values = [list(rows) for rows in array]

    def __getitem__(self, item):
        return self.values[item]

    def __add__(self, other):
        return Matrix(mat_sum(self.values, other))

    def __sub__(self, other):
        return Matrix(mat_sum(self.values, mat_scale(other, -1)))

    def __mul__(self, other):
        return Matrix(mat_product(self.values, other))

    def __str__(self):
        s = '\t['
        for r in self.values:
            s += str(r) + ',\n\t '
        return s[:-4] + ']'

    def __repr__(self):
        return str(self)

    def apply(self, f):
        v = self.values
        out = [[f(e) for e in r] for r in v]
        return Matrix(out)

    def transpose(self):
        return Matrix(transpose(self.values))

    def inverse(self):
        return Matrix(getMatrixInverse(self.values))

    def det(self):
        return getMatrixDeternminant(self.values)


class Quaternion:
    def __init__(self, w, x, y, z):
        """
        A Quaternion object
        :param w: q0
        :param x: q1
        :param y: q2
        :param z: q3
        """
        assert isinstance(w, (int, float)), "All components must be real numbers"
        assert isinstance(x, (int, float)), "All components must be real numbers"
        assert isinstance(y, (int, float)), "All components must be real numbers"
        assert isinstance(z, (int, float)), "All components must be real numbers"

        self.w = w
        self.x = x
        self.y = y
        self.z = z

    def components(self):
        return self.w, self.x, self.y, self.z

    def tb(self):
        """
        Wrapper method of the quaternion2tb method from the 'quaternions.py' file
        :return: A (roll, pitch, yaw) tuple corresponding to the Tait-Bryan angles - in radians -.
        """
        angle = quaternion2tb(self.w, self.x, self.y, self.z)
        return angle

    def aa(self):
        """
        Wrapper method of the quaternion2aa method from the 'quaternions.py' file
        :return: An (angle, axis) tuple - in radians -.
        """
        angle, axis = quaternion2aa(self.w, self.x, self.y, self.z)
        return angle, axis

    def norm(self):
        n = self.w ** 2
        n += self.x ** 2
        n += self.y ** 2
        n += self.z ** 2
        return n ** 0.5

    def normalize(self):
        n = self.norm()
        if n > 0:
            self.w /= n
            self.x /= n
            self.y /= n
            self.z /= n

    def conjugate(self):
        w = self.w
        x = -1 * self.x
        y = -1 * self.y
        z = -1 * self.z
        return Quaternion(w, x, y, z)

    def inverse(self):
        n = self.norm()
        w, x, y, z = vector_scale(self.conjugate().components(), n ** -2)
        return Quaternion(w, x, y, z)

    def dot_product(self, other):
        if not isinstance(other, (Quaternion, tuple)):
            raise TypeError("Cannot perform dot product between a Quaternion and a(n) {}".format(type(other)))
        if isinstance(other, Quaternion):
            other = other.components()
        return dot_product(self.components(), other)

    def rotate(self, other):
        if not isinstance(other, Quaternion):
            raise TypeError("Cannot rotate a Quaternion by a non-Quaternion object")
        return other * (self * other.conjugate())

    def __add__(self, other):
        if not isinstance(other, Quaternion):
            raise TypeError("Cannot add a Quaternion object with a(n) {}".format(type(other)))
        w = self.w + other.w
        x = self.x + other.x
        y = self.y + other.y
        z = self.z + other.z
        return Quaternion(w, x, y, z)

    def __sub__(self, other):
        if not isinstance(other, Quaternion):
            raise TypeError("Cannot subtract a(n) {} from a Quaternion object".format(type(other)))
        w = self.w - other.w
        x = self.x - other.x
        y = self.y - other.y
        z = self.z - other.z
        return Quaternion(w, x, y, z)

    def __mul__(self, other):
        if not isinstance(other, (Quaternion, int, float)):
            raise TypeError("Cannot multiply a Quaternion object with a(n) {}".format(type(other)))

        if isinstance(other, Quaternion):
            a = self.components()
            b = other.components()

            r1, r2 = a[0], b[0]
            v1, v2 = a[1:], b[1:]
            w = r1 * r2 - dot_product(v1, v2)
            xyz = vector_scale(r1, v2)
            xyz = vector_sum(xyz, vector_scale(r2, v1))
            x, y, z = vector_sum(xyz, cross_product(v1, v2))
        else:
            a = self.components()
            w, x, y, z = vector_scale(a, other)
        return Quaternion(w, x, y, z)

    def __str__(self):
        return str(self.components())

    def __repr__(self):
        return str(self)


def quaternion2tb(w, x, y, z):
    """
    This method takes the quaternion components and returns all Tait-Bryan angles.

    The coordinate system has:
        Roll axis pointing forward,
        Pitch axis pointing to the left and,
        Yaw axis pointing upward

    Taken from https://math.stackexchange.com/questions/2975109/how-to-convert-euler-angles-to-quaternions-and-get-the-same-euler-angles-back-fr)

    :param w: q0
    :param x: q1
    :param y: q2
    :param z: q3
    :return: A (yaw, pitch, roll) tuple corresponding to the Tait-Bryan angles in radians
    """
    t0 = +2.0 * (w * x + y * z)
    t1 = +1.0 - 2.0 * (x * x + y * y)
    roll = atan2(t0, t1)
    t2 = +2.0 * (w * y - z * x)
    t2 = +1.0 if t2 > +1.0 else t2
    t2 = -1.0 if t2 < -1.0 else t2
    pitch = asin(t2)
    t3 = +2.0 * (w * z + x * y)
    t4 = +1.0 - 2.0 * (y * y + z * z)
    yaw = atan2(t3, t4)
    return yaw, pitch, roll


def tb2quaternion(yaw, pitch, roll):
    """
    This method takes the Tait-Bryan angles (in radians) and returns a quaternion object.

    The coordinate system performs Yaw -> Pitch -> Roll to compute the
    resulting quaternion.

    Taken from https://math.stackexchange.com/questions/2975109/how-to-convert-euler-angles-to-quaternions-and-get-the-same-euler-angles-back-fr)

    :param yaw: Angle respective to the axis pointing upward (z-axis)
    :param pitch: Angle respective to the axis pointing to the left (y-axis)
    :param roll: Angle respective to the axis pointing forward (x-axis)
    :return: A Quaternion object
    """
    cy = cos(yaw * 0.5)
    sy = sin(yaw * 0.5)
    cp = cos(pitch * 0.5)
    sp = sin(pitch * 0.5)
    cr = cos(roll * 0.5)
    sr = sin(roll * 0.5)

    w = cy * cp * cr + sy * sp * sr
    x = cy * cp * sr - sy * sp * cr
    y = sy * cp * sr + cy * sp * cr
    z = sy * cp * cr - cy * sp * sr

    return Quaternion(w, x, y, z)


def aa2quaternion(angle, axis):
    """
    This method converts an angle-axis representation to a quaternion
    :param angle: Object angle respective to the axis (in radians)
    :param axis: Vector representing the axis in 3D space
    :return: A quaternion object
    """
    assert isinstance(angle, (int, float)), "Angle must be a real number (in radians)"
    assert len(axis) == 3, "Axis must be a 3D vector"

    w = cos(angle * 0.5)
    s = sin(angle * 0.5)
    x, y, z = vector_scale(axis, s)
    return Quaternion(w, x, y, z)


def quaternion2aa(w, x, y, z):
    """
    This method takes the quaternion components and return their angle-axis representation
    :param w: q0
    :param x: q1
    :param y: q2
    :param z: q3
    :return: An (angle, axis) tuple (angle in radians)
    """
    angle = 2 * acos(w)
    s = (1 - w ** 2) ** 0.5
    axis = vector_scale((x, y, z), s)
    return angle, axis


def vector_sum(a, b):
    if isinstance(a, (int, float)):
        if isinstance(b, (int, float)):
            return a + b
        else:
            assert isinstance(b, (list, tuple)), "Inputs must be one dimensional vectors"
            assert isinstance(b[0], (int, float)), "Inputs must be one dimensional vectors"
            return [eb + a for eb in b]
    else:
        assert isinstance(a, (list, tuple)), "Inputs must be one dimensional vectors"
        assert isinstance(a[0], (int, float)), "Inputs must be one dimensional vectors"

        if isinstance(b, (int, float)):
            return [ea + b for ea in a]
        else:
            assert isinstance(b, (list, tuple)), "Inputs must be one dimensional vectors"
            assert isinstance(b[0], (int, float)), "Inputs must be one dimensional vectors"
            assert len(a) == len(b), "The vectors to sum must be the same length"
            r = []
            for ea, eb in zip(a, b):
                r.append(ea+eb)
            return r


def vector_norm(a):
    n = dot_product(a, a)
    return n ** 0.5


def dot_product(a, b):
    assert len(a) == len(b), "The vectors to multiply must be the same length"
    p = 0
    for ea, eb in zip(a, b):
        p += ea * eb
    return p


def cross_product(a, b):
    assert len(a) == 3, "Both vectors must be 3D vectors"
    assert len(b) == 3, "Both vectors must be 3D vectors"
    c = (a[1]*b[2] - a[2]*b[1],
         a[2]*b[0] - a[0]*b[2],
         a[0]*b[1] - a[1]*b[0])
    return c


def vector_scale(a, b):
    assert isinstance(a, (int, float)) or isinstance(b, (int, float)), "One of the elements must be a real number"

    if isinstance(a, (int, float)):
        if isinstance(b, (int, float)):
            return a * b
        else:
            r = [a*eb for eb in b]
            return r
    else:
        r = [b * ea for ea in a]
        return r


def mat_scale(a, b):
    if isinstance(a, Matrix):
        a = a.values
    if isinstance(b, Matrix):
        b = b.values
    assert isinstance(a, (list, tuple)) and isinstance(a[0], (list, tuple)), "First input must be a matrix"
    assert isinstance(b, (int, float)), "Second input must be a real number"
    a = [vector_scale(ra, b) for ra in a]
    return a


def mat_sum(a, b):
    if isinstance(a, Matrix):
        a = a.values
    if isinstance(b, Matrix):
        b = b.values
    assert isinstance(a, (list, tuple)) and isinstance(a[0], (list, tuple)), "First input must be a matrix"
    assert isinstance(b, (list, tuple)) and isinstance(b[0], (list, tuple)), "Second input must be a matrix"
    assert len(a) == len(b) and len(a[0]) == len(b[0]), "Matrices should be the same size"
    out = [vector_sum(ea, eb) for ea, eb in zip(a, b)]
    return out


def mat_product(a, b):
    if isinstance(a, Matrix):
        a = a.values
    if isinstance(b, Matrix):
        b = b.values
    if isinstance(a, (int, float)) or isinstance(b, (int, float)):
        if isinstance(a[0], (int, float)):
            return vector_scale(a, b)
        else:
            return mat_scale(a, b)
    else:
        if isinstance(a[0], (int, float)):
            if isinstance(b[0], (int, float)):
                # 1D Row Vector X 1D Column vector
                return dot_product(a, b)
            else:
                # 1D Row Vector X 2D Matrix
                assert isinstance(b[0], (list, tuple)), "The 'b' Matrix should be a list-of-lists (or tuples)"
                b = transpose(b)
                return [dot_product(a, eb) for eb in b]
        else:
            # 2D Matrix X 2D Matrix
            assert isinstance(a[0], (list, tuple)), "The 'a' Matrix should be a list-of-lists (or tuples)"
            if isinstance(b[0], (int, float)):
                # 2D Matrix X 1D Column vector
                return [dot_product(ea, b) for ea in a]
            else:
                assert isinstance(b[0], (list, tuple)), "The 'b' Matrix should be a list-of-lists (or tuples)"
                assert len(b) == len(a[0]), "The number of columns in 'a' should be equal to the number of rows in 'b'"
                b = transpose(b)
                out = []
                for ra in a:
                    eo = [dot_product(ra, cb) for cb in b]
                    out.append(eo)
                return out


def zeros(a):
    assert isinstance(a, (int, float)), "The input must be a real number"
    m = []
    for i in range(a):
        r = [0 for j in range(a)]
        m.append(r)
    return Matrix(m)


def identity(a):
    assert isinstance(a, (int, float)), "The input must be a real number"
    m = zeros(a)
    for i in range(a):
        m[i][i] = 1
    return m


def transpose(a):
    return [list(i) for i in zip(*a)]


def hamilton_product(a, b):
    assert isinstance(a, (Quaternion, tuple)) and \
           isinstance(b, (Quaternion, tuple)), "Both elements must be either Quaternion or a tuple"
    if isinstance(a, Quaternion):
        a = a.components()
    else:
        for ea in a:
            assert isinstance(ea, (int, float)), "All tuple elements must be real numbers"
    if isinstance(b, Quaternion):
        b = b.components()
    else:
        for eb in b:
            assert isinstance(eb, (int, float)), "All tuple elements must be real numbers"
    assert len(a) == 4 and len(b) == 4, "Both objects must have 4 elements"

    a1, b1, c1, d1 = a
    a2, b2, c2, d2 = b

    w = a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2
    x = a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2
    y = a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2
    z = a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2

    return Quaternion(w, x, y, z)


# Third party code -----------------------------------------------------------------------------------------------------
# Taken from (https://stackoverflow.com/questions/32114054/matrix-inversion-without-numpy/39881366)

def getMatrixMinor(m,i,j):
    return [row[:j] + row[j+1:] for row in (m[:i]+m[i+1:])]


def getMatrixDeternminant(m):
    #  base case for 2x2 matrix
    if len(m) == 2:
        return m[0][0]*m[1][1]-m[0][1]*m[1][0]

    determinant = 0
    for c in range(len(m)):
        determinant += ((-1)**c)*m[0][c]*getMatrixDeternminant(getMatrixMinor(m,0,c))
    return determinant


def getMatrixInverse(m):
    determinant = getMatrixDeternminant(m)
    #  special case for 2x2 matrix:
    if len(m) == 2:
        return [[m[1][1]/determinant, -1*m[0][1]/determinant],
                [-1*m[1][0]/determinant, m[0][0]/determinant]]

    #  find matrix of cofactors
    cofactors = []
    for r in range(len(m)):
        cofactorRow = []
        for c in range(len(m)):
            minor = getMatrixMinor(m, r, c)
            cofactorRow.append(((-1)**(r+c)) * getMatrixDeternminant(minor))
        cofactors.append(cofactorRow)
    cofactors = transpose(cofactors)
    for r in range(len(cofactors)):
        for c in range(len(cofactors)):
            cofactors[r][c] = cofactors[r][c]/determinant
    return cofactors
