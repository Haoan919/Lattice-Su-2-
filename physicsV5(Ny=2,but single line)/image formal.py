import numpy as np
from numba import njit
import matplotlib.pyplot as plt

from functions import *

sx = np.array(((0, 1), (1, 0)), complex)
sy = np.array(((0, -1j), (1j, 0)), complex)
sz = np.array(((1, 0), (0, -1)), complex)

su2 = 2

Ny = 2


def initialize_fields(start, Nt, Nx):
    """Initialize gauge configuration and Higgs field"""

    U = np.zeros((Nt, Nx, Ny, 3, su2, su2), complex)

    for t in range(Nt):
        for x in range(Nx):
            for y in range(Ny):
                for mu in range(3):
                    if y == 1 and mu == 2:
                        U[t, x, y, mu] = U[t, x, y - 1, mu].conj().T
                    else:
                        if start == 0:
                            U[t, x, y, mu] = np.identity(su2)
                        if start == 1:
                            U[t, x, y, mu] = SU2SingleMatrix()
    return U



@njit()
def HeatBath_updating_links(U, beta, Nt, Nx):
    for t in range(Nt):
        for x in range(Nx):
            for y in range(Ny):
                for mu in range(3):
                    staple = staple_calculus(t, x, y, mu, U, Nt, Nx)
                    U[t, x, y, mu] = HB_gauge(staple, beta)

    return U


################################ UPDATING ALGORITHMS ################################


@njit()
def Metropolis(U, beta, hits, Nt, Nx):
    """Execution of Metropolis, checking every single site 10 times."""

    for t in range(Nt):
        for x in range(Nx):
            for y in range(Ny):
                for mu in range(3):

                    staple = staple_calculus(t, x, y, mu, U, Nt, Nx)

                    for _ in range(hits):

                        old_link = U[t, x, y, mu].copy()
                        S_old = calculate_S(old_link, staple, beta)

                        su2matrix = SU2SingleMatrix()
                        new_link = np.dot(su2matrix, old_link)
                        S_new = calculate_S(new_link, staple, beta)
                        dS = S_new - S_old

                        if dS < 0:
                            U[t, x, y, mu] = new_link
                        else:
                            if np.exp(-dS) > np.random.uniform(0, 1):
                                U[t, x, y, mu] = new_link
                            else:
                                U[t, x, y, mu] = old_link
    return U


@njit()
def OverRelaxation(U, Nt, Nx):
    for t in range(Nt):
        for x in range(Nx):
            for y in range(Ny):
                for mu in range(3):

                    A = staple_calculus(t, x, y, mu, U, Nt, Nx)

                    a = np.sqrt((np.linalg.det(A)))

                    Utemp = U[t, x, y, mu]

                    if a.real != 0:

                        V = A / a
                        Uprime = V.conj().T @ Utemp.conj().T @ V.conj().T
                        # print(np.linalg.det(Uprime))
                        U[t, x, y, mu] = Uprime

                    else:

                        U[t, x, y, mu] = SU2SingleMatrix()

    return U


##########################################################################################à

if __name__ == "__main__":

    import time

    start = time.time()

    Nt = 20
    Nx = 10
    measures = 20
    beta_vec = np.linspace(0.1, 10, 20).tolist()

    U = initialize_fields(1, Nt, Nx)

    obs11 = []
    obs22 = []
    obs33 = []

    heatbath = True
    metropolis = False
    overrelax = True

    for b in beta_vec:
        print("exe for beta ", b)

        smean11 = []
        smean22 = []
        smean33 = []

        for m in range(measures):

            if heatbath:
                U = HeatBath_updating_links(U, b, Nt, Nx)

            if metropolis:
                U = Metropolis(U, b, hits=10)

            if overrelax:

                for _ in range(2):
                    U = OverRelaxation(U, Nt, Nx)

            temp11 = WilsonAction(0, 1, 1, 1, U, Nt, Nx)
            temp22 = WilsonAction(0, 1, 2, 2, U, Nt, Nx)
            temp33 = WilsonAction(0, 1, 3, 3, U, Nt, Nx)
            print(temp11)

            smean11.append(temp11)
            smean22.append(temp22)
            smean33.append(temp33)

        obs11.append(np.mean(smean11))
        obs22.append(np.mean(smean22))
        obs33.append(np.mean(smean33))

    print("tempo:", round(time.time() - start, 2))

    plt.figure()
    plt.plot(beta_vec, -np.log(obs11), 'go')
    plt.plot(beta_vec, -np.log(obs22) / 4, 'ro')
    plt.plot(beta_vec, -np.log(obs33) / 9, 'bo')
    plt.show()

    plt.figure()
    plt.plot(beta_vec, obs11,'g.')
    plt.plot(beta_vec, obs22,'b.')
    plt.plot(beta_vec, obs33,'b.')
    plt.show()