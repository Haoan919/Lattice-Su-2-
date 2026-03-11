import numpy as np
from numba import njit

sx = np.array(((0, 1), (1, 0)), complex)
sy = np.array(((0, -1j), (1j, 0)), complex)
sz = np.array(((1, 0), (0, -1)), complex)

su2 = 2


@njit()
def WilsonAction12(R, T, U, Nt, Nx, Ny):
    somma = 0
    for t in range(Nt):
        for x in range(0, 2 * Nx, 2):
            for y in range(Ny):
                loop = np.identity(su2) + 0j

                for i in range(2 * R):
                    loop = np.dot(loop, U[t, (x + i) % (2 * Nx), y, 1])

                for j in range(T - 1):
                    loop = np.dot(loop, U[t, (x + 2 * R) % (2 * Nx), (y + j) % Ny, 2])
                    loop = np.dot(loop, U[t, (x + 2 * R - 1) % (2 * Nx), (y + j + 1) % Ny, 1])

                loop = np.dot(loop, U[t, (x + 2 * R) % (2 * Nx), (y + T - 1) % Ny, 2])

                for i in range(2 * R - 2, -2, -1):
                    loop = np.dot(loop, U[t, (x + i) % (2 * Nx), (y + T) % Ny, 1].conj().T)

                for j in range(T - 1, 0, -1):
                    loop = np.dot(loop, U[t, x, (y + j) % Ny, 2].conj().T)
                    loop = np.dot(loop, U[t, (x - 1) % (2 * Nx), (y + j) % Ny, 1].conj().T)
                loop = np.dot(loop, U[t, x, y, 2].conj().T)

                somma += np.trace(loop).real / su2

    return somma / (Nt * Nx * Ny)


@njit()
def WilsonAction01(R, T, U, Nt, Nx, Ny):
    somma = 0
    for t in range(Nt):
        for x in range(0, 2 * Nx, 2):
            for y in range(Ny):
                loop = np.identity(su2) + 0j

                for i in range(R):
                    loop = np.dot(loop, U[(t + i) % Nt, x, y, 0])

                for j in range(2 * T):
                    loop = np.dot(loop, U[(t + R) % Nt, (x + i) % (2 * Nx), y, 1])

                for i in range(R - 1, -1, -1):
                    loop = np.dot(loop, U[(t + i) % Nt, (x + 2 * T) % (2 * Nx), y, 0].conj().T)

                for j in range(2 * T - 1, -1, -1):
                    loop = np.dot(loop, U[t, (x + i) % (2 * Nx), y, 1].conj().T)

                somma += np.trace(loop).real / su2

    return somma / (Nt * Nx * Ny)


@njit()
def WilsonAction02(R, T, U, Nt, Nx, Ny):
    somma = 0
    for t in range(Nt):
        for x in range(0, 2 * Nx, 2):
            for y in range(Ny):
                loop = np.identity(su2) + 0j

                for i in range(R):
                    loop = np.dot(loop, U[(t + i) % Nt, x, y, 0])

                for j in range(T - 1):
                    loop = np.dot(loop, U[(t + R) % Nt, x, (y + j) % Ny, 2])
                    loop = np.dot(loop, U[(t + R) % Nt, (x - 1) % (2 * Nx), (y + j + 1) % Ny, 1])
                loop = np.dot(loop, U[(t + R) % Nt, x, (y + T - 1) % Ny, 2])

                for i in range(R - 1, -1, -1):
                    loop = np.dot(loop, U[(t + i) % Nt, (x - 1) % (2 * Nx), (y + T) % Ny, 0].conj().T)

                for j in range(T - 1, 0, -1):
                    loop = np.dot(loop, U[t, x, (y + j) % Ny, 2].conj().T)
                    loop = np.dot(loop, U[t, (x - 1) % (2 * Nx), (y + j) % Ny, 1].conj().T)
                loop = np.dot(loop, U[t, x, y, 2].conj().T)

                somma += np.trace(loop).real / su2

    return somma / (Nt * Nx * Ny)


@njit()
def staple_calculus(t, x, y, mu, U, Nt, Nx, Ny):
    staple_start = np.array(((0 + 0j, 0 + 0j), (0 + 0j, 0 + 0j)))

    if mu == 0:
        staple_start += (U[(t + 1) % Nt, x, y, 1] @
                         U[t, (x + 1) % (2 * Nx), y, 0].conj().T @
                         U[t, x, y, 1].conj().T)
        staple_start += (U[(t + 1) % Nt, (x - 1) % (2 * Nx), y, 1].conj().T @
                         U[t, (x - 1) % (2 * Nx), y, 0].conj().T @
                         U[t, (x - 1) % (2 * Nx), y, 1])
        if x % 2 == 0:
            staple_start += (U[(t + 1) % Nt, x, y, 2] @
                             U[t, (x - 1) % (2 * Nx), (y + 1) % Ny, 0].conj().T @
                             U[t, x, y, 2].conj().T)
        else:
            staple_start += (U[(t + 1) % Nt, (x + 1) % (2 * Nx), (y - 1) % Ny, 2].conj().T @
                             U[t, (x + 1) % (2 * Nx), (y - 1) % Ny, 0].conj().T @
                             U[t, (x + 1) % (2 * Nx), (y - 1) % Ny, 2])
    if mu == 2:
        staple_start += (U[t, (x - 1) % (2 * Nx), (y + 1) % Ny, 0] @
                         U[(t + 1) % Nt, x, y, 2].conj().T @
                         U[t, x, y, 0].conj().T)
        staple_start += (U[(t - 1) % Nt, (x - 1) % (2 * Nx), (y + 1) % Ny, 0].conj().T @
                         U[(t - 1) % Nt, x, y, 2].conj().T @
                         U[(t - 1) % Nt, x, y, 0])
        staple_start += (U[t, (x - 1) % (2 * Nx), (y + 1) % Ny, 1] @
                         U[t, x, (y + 1) % Ny, 1] @
                         U[t, (x + 2) % (2 * Nx), y, 2].conj().T @
                         U[t, (x + 1) % (2 * Nx), y, 1].conj().T @
                         U[t, x, y, 1].conj().T)
        staple_start += (U[t, (x - 2) % (2 * Nx), (y + 1) % Ny, 1].conj().T @
                         U[t, (x - 3) % (2 * Nx), (y + 1) % Ny, 1].conj().T @
                         U[t, (x - 2) % (2 * Nx), y, 2].conj().T @
                         U[t, (x - 2) % (2 * Nx), y, 1] @
                         U[t, (x - 1) % (2 * Nx), y, 1])

    if mu == 1:
        staple_start += (U[t, (x + 1) % (2 * Nx), y, 0] @
                         U[(t + 1) % Nt, x, y, 1].conj().T @
                         U[t, x, y, 0].conj().T)
        staple_start += (U[(t - 1) % Nt, (x + 1) % (2 * Nx), y, 0].conj().T @
                         U[(t - 1) % Nt, x, y, 1].conj().T @
                         U[(t - 1) % Nt, x, y, 0])
        if x % 2 == 0:
            staple_start += (U[t, (x + 1) % (2 * Nx), y, 1] @
                             U[t, (x + 2) % (2 * Nx), y, 2] @
                             U[t, x, (y + 1) % Ny, 1].conj().T @
                             U[t, (x - 1) % (2 * Nx), (y + 1) % Ny, 1].conj().T @
                             U[t, x, y, 2].conj().T)
            staple_start += (U[t, (x + 2) % (2 * Nx), (y - 1) % Ny, 2].conj().T @
                             U[t, (x + 1) % (2 * Nx), (y - 1) % Ny, 1].conj().T @
                             U[t, x, (y - 1) % Ny, 1].conj().T @
                             U[t, x, (y - 1) % Ny, 2] @
                             U[t, (x - 1) % (2 * Nx), y, 1])
        else:
            staple_start += (U[t, (x + 1) % (2 * Nx), y, 2] @
                             U[t, (x - 1) % (2 * Nx), (y + 1) % Ny, 1].conj().T @
                             U[t, (x - 2) % (2 * Nx), (y + 1) % Ny, 1].conj().T @
                             U[t, (x - 1) % (2 * Nx), y, 2].conj().T @
                             U[t, (x - 1) % (2 * Nx), y, 1])
            staple_start += (U[t, (x + 1) % (2 * Nx), y, 1] @
                             U[t, (x + 3) % (2 * Nx), (y - 1) % Ny, 2].conj().T @
                             U[t, (x + 2) % (2 * Nx), (y - 1) % Ny, 1].conj().T @
                             U[t, (x + 1) % (2 * Nx), (y - 1) % Ny, 1].conj().T @
                             U[t, (x + 1) % (2 * Nx), (y - 1) % Ny, 2])

    return staple_start


@njit()
def SU2SingleMatrix(epsilon=0.2):
    r0 = np.random.uniform(-0.5, 0.5)
    x0 = np.sign(r0) * np.sqrt(1 - epsilon ** 2)

    r = np.random.random((3)) - 0.5
    x = epsilon * r / np.linalg.norm(r)

    SU2Matrix = x0 * np.identity(2) + 1j * x[0] * sx + 1j * x[1] * sy + 1j * x[2] * sz

    return SU2Matrix


@njit()
def HB_gauge(staple, beta, kind=2):
    if kind == 1:

        detStaple = np.linalg.det(staple)
        alpha = np.sqrt(detStaple)

        r1 = np.random.uniform(0, 1)
        r2 = np.random.uniform(0, 1)
        r3 = np.random.uniform(0, 1)

        lambda2 = (
                -1
                / (2 * alpha * beta)
                * (np.log(r1) + np.cos(2 * np.pi * r2) ** 2 * np.log(r3))
        )

        while (np.random.uniform(0, 1)) ** 2 > (1 - lambda2.real):
            r1 = np.random.uniform(0, 1)
            r2 = 1 - np.random.uniform(0, 1)
            r3 = 1 - np.random.uniform(0, 1)

            lambda2 = (
                    -1
                    / (2 * alpha * beta)
                    * (np.log(r1) + np.cos(2 * np.pi * r2) ** 2 * np.log(r3))
            )

        x0 = 1 - 2 * lambda2

        a0 = x0

        a1 = np.random.uniform(-1, 1)
        a2 = np.random.uniform(-1, 1)
        a3 = np.random.uniform(-1, 1)

        while np.sqrt((a1 ** 2 + a2 ** 2 + a3 ** 2)) > 1:
            a1 = np.random.uniform(-1, 1)
            a2 = np.random.uniform(-1, 1)
            a3 = np.random.uniform(-1, 1)

        r = np.sqrt(1 - a0 ** 2)
        norm = np.sqrt(a1 ** 2 + a2 ** 2 + a3 ** 2)
        a1 = a1 * r / norm
        a2 = a2 * r / norm
        a3 = a3 * r / norm
        # Ulink = np.array(((a0 + 1j * a3, a2 + 1j * a1), (-a2 + 1j * a1, a0 - 1j * a3)))
        Ulink = a0 * np.identity(su2) + 1j * a1 * sx + 1j * a2 * sy + 1j * a3 * sz

        return Ulink

    if kind == 2:

        w = normalize(getA(staple))
        a = np.sqrt(np.abs(np.linalg.det(staple)))
        wbar = quaternion((w))

        if a != 0:
            xw = quaternion(sampleA(beta, a))

            xx = xw @ wbar.conj().T  ###!!!!warning!!!

            return xx

        else:
            return SU2SingleMatrix()


@njit()
def getA(W):
    """Construct the vector needed for the quaternion"""

    a0 = ((W[0, 0] + W[1, 1])).real / 2
    a1 = ((W[0, 1] + W[1, 0])).imag / 2
    a2 = ((W[0, 1] - W[1, 0])).real / 2
    a3 = ((W[0, 0] - W[1, 1])).imag / 2
    # Avector = [a0, a1, a2, a3]
    Avector = np.array((a0, a1, a2, a3))

    return Avector


@njit()
def sampleA(a, beta):
    """choose a0 with P(a0) ~ sqrt(1 - a0^2) * exp(beta * k * a0)"""

    w = np.exp(-2 * beta * a)
    xtrial = np.random.uniform(0, 1) * (1 - w) + w
    a0 = 1 + np.log(xtrial) / (beta * a)

    while (1 - a0 ** 2) < np.random.uniform(0, 1):
        xtrial = np.random.uniform(0, 1) * (1 - w) + w
        a0 = 1 + np.log(xtrial) / (beta * a)

    r = np.sqrt(1 - a0 ** 2)
    a1 = np.random.normal()
    a2 = np.random.normal()
    a3 = np.random.normal()

    while (a1 ** 2 + a2 ** 2 + a3 ** 2) > 1:
        a1 = np.random.normal()
        a2 = np.random.normal()
        a3 = np.random.normal()

    norm = np.sqrt(a1 ** 2 + a2 ** 2 + a3 ** 2)

    a1 = a1 * r / norm
    a2 = a2 * r / norm
    a3 = a3 * r / norm

    avec = np.array((a0, a1, a2, a3))
    avec = normalize(avec)

    return avec


@njit()
def quaternion(vec):
    """produces quaternion from a vector of complex and real numbers"""

    a11 = vec[0] + vec[3] * 1j
    a12 = vec[2] + vec[1] * 1j
    a21 = -vec[2] + vec[1] * 1j
    a22 = vec[0] - vec[3] * 1j

    quat = np.array(((a11, a12), (a21, a22)))

    return quat


@njit()
def calculate_S(link, stapla, beta):
    ans = np.trace(np.dot(stapla, link)).real
    return -beta * ans / su2


@njit()
def normalize(v):
    lun = len(v)
    s = 0

    for i in range(lun):
        s += v[i] ** 2

    return v / np.sqrt(s)


@njit()
def GramSchmidt(A, exe):
    """GS orthogonalization, if exe==False => the function normalizes only;
    non dovrebbe servire a niente"""

    n = len(A)

    if exe:
        A[:, 0] = normalize(A[:, 0])

        for i in range(1, n):
            Ai = A[:, i]
            for j in range(0, i):
                Aj = A[:, j]
                t = Ai.dot(Aj)
                Ai = Ai - t * Aj
            A[:, i] = normalize(Ai)
    else:
        for k in range(n):
            # print(np.linalg.norm(A[:, k]))
            A[:, k] = A[:, k] / np.linalg.norm(A[:, k])

    return A