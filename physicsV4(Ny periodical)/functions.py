import numpy as np
from numba import njit

sx = np.array(((0, 1), (1, 0)), complex)
sy = np.array(((0, -1j), (1j, 0)), complex)
sz = np.array(((1, 0), (0, -1)), complex)

su2 = 2

Ny = 1


@njit()
def WilsonAction(nu, mu, R, T, U, Nt, Nx):

    """Name says"""

    somma = 0
    for t in range(Nt):
        for x in range(Nx):
            for y in range(Ny):

                a_nu = [0, 0, 0]
                a_nu[nu] = 1
                a_mu = [0, 0, 0]
                a_mu[mu] = 1
                i = 0
                j = 0

                loop = np.identity(su2) + 0j

                a_nu = [0, 0, 0]
                a_nu[nu] = 1

                for i in range(R):

                    loop = np.dot(
                        loop,
                        U[
                            (t + i * a_mu[0]) % Nt,
                            (x + i * a_mu[1]) % Nx,
                            (y + i * a_mu[2]) % Ny,
                            mu,
                        ],
                    )

                for j in range(T):
                    loop = np.dot(
                        loop,
                        U[
                            (t + T * a_mu[0] + j * a_nu[0]) % Nt,
                            (x + T * a_mu[1] + j * a_nu[1]) % Nx,
                            (y + T * a_mu[2] + j * a_nu[2]) % Ny,
                            nu,
                        ],
                    )

                for i in range(R - 1, -1, -1):
                    loop = np.dot(
                        loop,
                        U[
                            (t + i * a_mu[0] + R * a_nu[0]) % Nt,
                            (x + i * a_mu[1] + R * a_nu[1]) % Nx,
                            (y + i * a_mu[2] + R * a_nu[2]) % Ny,
                            mu,
                        ]
                        .conj()
                        .T,
                    )

                for j in range(T - 1, -1, -1):
                    loop = np.dot(
                        loop,
                        U[
                            (t + j * a_nu[0]) % Nt,
                            (x + j * a_nu[1]) % Nx,
                            (y + j * a_nu[2]) % Ny,
                            nu,
                        ]
                        .conj()
                        .T,
                    )

                somma += np.trace(loop).real / su2

    return somma / (Nt * Nx * Ny)
@njit()
def staple_calculus(t, x, y, mu, U, Nt, Nx):

    """Calculate the contribution (interaction) of the three links sorrounding the link that we want to update"""
    # staple_start = np.zeros((su3, su3), complex)
    staple_start = np.array(((0 + 0j, 0 + 0j), (0 + 0j, 0 + 0j)))
    # staple_start = np.empty((3, 3)) + 0j

    a_mu = [0, 0, 0]
    a_mu[mu] = 1

    for nu in range(3):

        if mu != nu:
            # nu = 0
            # while nu < mu:
            a_nu = [0, 0, 0]
            a_nu[nu] = 1

            # product of 6 matrices
            staple_start += (
                U[
                    (t + a_mu[0]) % Nt,
                    (x + a_mu[1]) % Nx,
                    (y + a_mu[2]) % Ny,
                    nu,
                ]
                @ U[
                    (t + a_nu[0]) % Nt,
                    (x + a_nu[1]) % Nx,
                    (y + a_nu[2]) % Ny,
                    mu,
                ]
                .conj()
                .T
                @ U[t, x, y, nu].conj().T
            )

            staple_start += (
                U[
                    (t + a_mu[0] - a_nu[0]) % Nt,
                    (x + a_mu[1] - a_nu[1]) % Nx,
                    (y + a_mu[2] - a_nu[2]) % Ny,
                    nu,
                ]
                .conj()
                .T
                @ U[
                    (t - a_nu[0]) % Nt,
                    (x - a_nu[1]) % Nx,
                    (y - a_nu[2]) % Ny,
                    mu,
                ]
                .conj()
                .T
                @ U[
                    (t - a_nu[0]) % Nt,
                    (x - a_nu[1]) % Nx,
                    (y - a_nu[2]) % Ny,
                    nu,
                ]
            )
            # nu += 1
        else:
            continue

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