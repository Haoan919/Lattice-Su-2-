import numpy as np
from numba import njit
import matplotlib.pyplot as plt

from functions import *

sx = np.array(((0, 1), (1, 0)), complex)
sy = np.array(((0, -1j), (1j, 0)), complex)
sz = np.array(((1, 0), (0, -1)), complex)

su2 = 2
N = 5


def initialize_fields(start):

    """Initialize gauge configuration and Higgs field"""

    U = np.zeros((N, N, N, N, 4, su2, su2), complex)
    UHiggs = U.copy()

    for t in range(N):
        for x in range(N):
            for y in range(N):
                for z in range(N):
                    for mu in range(4):

                        UHiggs[t, x, y, z, mu] = SU2SingleMatrix()

                        if start == 0:
                            U[t, x, y, z, mu] = np.identity(su2)
                        if start == 1:
                            U[t, x, y, z, mu] = SU2SingleMatrix()

    return U, UHiggs


@njit()
def HeatBath_updating_links(U, beta):

    for t in range(N):
        for x in range(N):
            for y in range(N):
                for z in range(N):
                    for mu in range(4):
                        staple = staple_calculus(t, x, y, z, mu, U)
                        U[t, x, y, z, mu] = HB_gauge(staple, beta)

    return U


################################ UPDATING ALGORITHMS ################################


@njit()
def Metropolis(U, beta, hits):

    """Execution of Metropolis, checking every single site 10 times."""

    for t in range(N):
        for x in range(N):
            for y in range(N):
                for z in range(N):
                    for mu in range(4):

                        staple = staple_calculus(t, x, y, z, mu, U)

                        for _ in range(hits):

                            old_link = U[t, x, y, z, mu].copy()
                            S_old = calculate_S(old_link, staple, beta)

                            su2matrix = SU2SingleMatrix()
                            new_link = np.dot(su2matrix, old_link)
                            S_new = calculate_S(new_link, staple, beta)
                            dS = S_new - S_old

                            if dS < 0:
                                U[t, x, y, z, mu] = new_link
                            else:
                                if np.exp(-dS) > np.random.uniform(0, 1):
                                    U[t, x, y, z, mu] = new_link
                                else:
                                    U[t, x, y, z, mu] = old_link
    return U


@njit()
def OverRelaxation(U):

    for t in range(N):
        for x in range(N):
            for y in range(N):
                for z in range(N):
                    for mu in range(4):

                        A = staple_calculus(t, x, y, z, mu, U)

                        a = np.sqrt((np.linalg.det(A)))

                        Utemp = U[t, x, y, z, mu]

                        if a.real != 0:

                            V = A / a
                            Uprime = V.conj().T @ Utemp.conj().T @ V.conj().T
                            # print(np.linalg.det(Uprime))
                            U[t, x, y, z, mu] = Uprime

                        else:

                            U[t, x, y, z, mu] = SU2SingleMatrix()

    return U


##########################################################################################à

if __name__ == "__main__":

    import time

    start = time.time()
    U, Uh = initialize_fields(1)
    measures = 40

    beta_vec = np.linspace(0.1, 10, 100).tolist()

    U, _ = initialize_fields(1)

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

                U = HeatBath_updating_links(U, b)

            if metropolis:

                U = Metropolis(U, b, hits=10)

            if overrelax:

                for _ in range(2):

                    U = OverRelaxation(U)

            temp11 = WilsonAction(1, 1, U)
            temp22 = WilsonAction(2, 2, U)
            temp33 = WilsonAction(3, 3, U)
            #wrong original write temp22 as 33
            print(temp11)

            smean11.append(temp11)
            smean22.append(temp22)
            smean33.append(temp33)

        obs11.append(np.mean(smean11))
        obs22.append(np.mean(smean22))
        obs33.append(np.mean(smean33))

    print("tempo:", round(time.time() - start, 2))

    plt.figure()
    plt.plot(beta_vec, obs11,'g.')
    plt.plot(beta_vec, obs22,'b.')
    plt.plot(beta_vec, obs33,'b.')
    plt.show()