import numpy as np
from numba import njit
import matplotlib.pyplot as plt

from functions import *

sx = np.array(((0, 1), (1, 0)), complex)
sy = np.array(((0, -1j), (1j, 0)), complex)
sz = np.array(((1, 0), (0, -1)), complex)

su2 = 2

Nx = 20
Ny = 2


def initialize_fields(start, Nt):
    """Initialize gauge configuration and Higgs field"""

    U = np.zeros((Nt, Nx, Ny, 3, su2, su2), complex)
    UHiggs = U.copy()

    for t in range(Nt):
        for x in range(Nx):
            for y in range(Ny):
                for mu in range(3):

                    UHiggs[t, x, y, mu] = SU2SingleMatrix()

                    if start == 0:
                        U[t, x, y, mu] = np.identity(su2)
                    if start == 1:
                        U[t, x, y, mu] = SU2SingleMatrix()

    return U, UHiggs


@njit()
def HeatBath_updating_links(U, beta, Nt):
    for t in range(Nt):
        for x in range(Nx):
            for y in range(Ny):
                for mu in range(3):
                    staple = staple_calculus(t, x, y, mu, U, Nt)
                    U[t, x, y, mu] = HB_gauge(staple, beta)

    return U


################################ UPDATING ALGORITHMS ################################


@njit()
def Metropolis(U, beta, hits, Nt):
    """Execution of Metropolis, checking every single site 10 times."""

    for t in range(Nt):
        for x in range(Nx):
            for y in range(Ny):
                for mu in range(3):

                    staple = staple_calculus(t, x, y, mu, U, Nt)

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
def OverRelaxation(U, Nt):
    for t in range(Nt):
        for x in range(Nx):
            for y in range(Ny):
                for mu in range(3):

                    A = staple_calculus(t, x, y, mu, U, Nt)

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
#这个下面写的是逐步递增Nt的情况

if __name__ == "__main__":

    import time

    start = time.time()

    measures = 40

    # beta_vec = np.linspace(0.1, 10, 20).tolist()

    b = 2.5
    Nt_vec = range(20, 40, 5)

    # U, _ = initialize_fields(1, Nt)

    obs11 = []

    heatbath = True
    metropolis = False
    overrelax = True

    for Nt in Nt_vec:
        print("exe for Nt ", Nt)

        U, _ = initialize_fields(1, Nt)

        smean11 = []

        for m in range(measures):

            if heatbath:
                U = HeatBath_updating_links(U, b, Nt)

            if metropolis:
                U = Metropolis(U, b, 10, Nt)

            if overrelax:

                for _ in range(2):
                    U = OverRelaxation(U, Nt)

            temp11 = WilsonAction(1, 1, U, Nt)
            print(temp11)

            smean11.append(temp11)

        obs11.append(np.mean(smean11))

    print("tempo:", round(time.time() - start, 2))

    plt.figure()
    plt.plot(Nt_vec, obs11, 'go-')
    plt.show()
#下面的这个讲的是不同wt对贝塔的函数，是for循环递增
    Nt = 20
    U, _ = initialize_fields(1, Nt)
    Wt_vec = range(1, 5, 1)
    obs1n = []

    for Wt in Wt_vec:
        print("exe for Wt ", Wt)

        smean11 = []

        for m in range(measures):

            if heatbath:
                U = HeatBath_updating_links(U, b, Nt)

            if metropolis:
                U = Metropolis(U, b, 10, Nt)

            if overrelax:

                for _ in range(2):
                    U = OverRelaxation(U, Nt)

            temp11 = WilsonAction(1, Wt, U, Nt)
            print(temp11)

            smean11.append(temp11)

        obs1n.append(np.mean(smean11))

    print("tempo:", round(time.time() - start, 2))

    plt.figure()
    plt.plot(Wt_vec, obs1n, 'go-')
    plt.show()