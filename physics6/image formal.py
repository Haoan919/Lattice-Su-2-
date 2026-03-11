import numpy as np
from numba import njit
import matplotlib.pyplot as plt

from functions import *

sx = np.array(((0, 1), (1, 0)), complex)
sy = np.array(((0, -1j), (1j, 0)), complex)
sz = np.array(((1, 0), (0, -1)), complex)

su2 = 2


def initialize_fields(start, Nt, Nx, Ny):
    """Initialize gauge configuration and Higgs field"""

    U = np.zeros((Nt, 2 * Nx, Ny, 3, su2, su2), complex)

    for t in range(Nt):
        for x in range(2 * Nx):
            for y in range(Ny):
                for mu in range(3):
                    if start == 0:
                        U[t, x, y, mu] = np.identity(su2)
                    if start == 1:
                        U[t, x, y, mu] = SU2SingleMatrix()
    return U


@njit()
def HeatBath_updating_links(U, beta, Nt, Nx, Ny):
    for t in range(Nt):
        for x in range(2 * Nx):
            for y in range(Ny):
                for mu in range(3):
                    staple = staple_calculus(t, x, y, mu, U, Nt, Nx, Ny)
                    U[t, x, y, mu] = HB_gauge(staple, beta)

    return U


################################ UPDATING ALGORITHMS ################################


@njit()
def Metropolis(U, beta, hits, Nt, Nx, Ny):
    """Execution of Metropolis, checking every single site 10 times."""

    for t in range(Nt):
        for x in range(2 * Nx):
            for y in range(Ny):
                for mu in range(3):

                    staple = staple_calculus(t, x, y, mu, U, Nt, Nx, Ny)

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
def OverRelaxation(U, Nt, Nx, Ny):
    for t in range(Nt):
        for x in range(2 * Nx):
            for y in range(Ny):
                for mu in range(3):

                    A = staple_calculus(t, x, y, mu, U, Nt, Nx, Ny)

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

    Nt = 5
    Nx = 3
    Ny = 3
    measures = 10
    beta_vec = np.linspace(0.1, 10, 20).tolist()

    U = initialize_fields(1, Nt, Nx, Ny)

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
                U = HeatBath_updating_links(U, b, Nt, Nx, Ny)

            if metropolis:
                U = Metropolis(U, b, hits=10)

            if overrelax:

                for _ in range(2):
                    U = OverRelaxation(U, Nt, Nx, Ny)

            temp11 = WilsonAction12(1, 1, U, Nt, Nx, Ny)
            temp22 = WilsonAction12(2, 2, U, Nt, Nx, Ny)
            temp33 = WilsonAction12(3, 3, U, Nt, Nx, Ny)
            print(temp11)

            smean11.append(temp11)
            smean22.append(temp22)
            smean33.append(temp33)

        obs11.append(np.mean(smean11))
        obs22.append(np.mean(smean22))
        obs33.append(np.mean(smean33))

    print("tempo:", round(time.time() - start, 2))

    plt.figure()
    plt.plot(beta_vec, -np.log(obs11),'go')
    plt.plot(beta_vec, -np.log(obs22)/4,'ro')
    plt.plot(beta_vec, -np.log(obs33)/9,'bo')
    plt.show()

    '''plt.figure()
    plt.plot(beta_vec, obs11, 'go')
    plt.plot(beta_vec, obs22, 'ro')
    plt.plot(beta_vec, obs33, 'bo')
    plt.show()'''
