from srrp import SRHDRiemannSolver
import matplotlib.pyplot as plt
import numpy as np


def simple_eval(rhoL, uxL, utL, pL, rhoR, uxR, utR, pR, gamma, T, x0=0.5):
    solver = SRHDRiemannSolver(rhoL, uxL, utL, pL, rhoR, uxR, utR, pR, gamma, x0, verbose=True)
    return solver


def complete_eval(rhoL, uxL, utL, pL, rhoR, uxR, utR, pR, gamma, T, x0=0.5, dx=0.5, N=1000):
    solver = SRHDRiemannSolver(rhoL, uxL, utL, pL, rhoR, uxR, utR, pR, gamma, x0, verbose=True)
    xs = np.linspace(x0 - dx, x0 + dx, N)

    qs = np.array([solver(x, T) for x in xs])

    regIdxs = np.array([solver.getRegIdx(x, T) for x in xs])

    for i in range(1,7):
        sel = regIdxs == i
        for j in range(4):
            plt.figure(j)
            plt.plot(xs[sel], qs[sel][:, j])

    for j in range(4):
        plt.figure(j)
        plt.xlim(np.min(xs), np.max(xs))
        y_min, y_max = np.min(qs[:, j]), np.max(qs[:, j])
        dy = y_max - y_min
        plt.ylim(y_min - 0.1*dy, y_max + 0.1*dy)

    W = qs[:,1]**2 + qs[:, 2]**2
    W = np.sqrt(1./(1. - W))

    #h = 1 + qs[:,3]/qs[:,0] * gamma/(gamma -1)
    #A = h * W * qs[:, 2]
    A = W


    plt.figure(4)
    for i in range(1,7):
        sel = regIdxs == i
        plt.plot(xs[sel], A[sel])

    plt.xlim(np.min(xs), np.max(xs))
    y_min, y_max = np.min(A), np.max(A)
    dy = y_max - y_min
    plt.ylim(y_min - 0.1 * dy, y_max + 0.1 * dy)

    return solver


def rezzolla_problem(n, T=0.4, invert=False, eval=simple_eval):
    print('rezzolla ', n)
    gamma = 5. / 3

    rho1 = 1
    p1 = 1
    p2 = 0.1
    rho2 = 0.125

    if n == 1:
        ux1, ux2 = 0.500, 0.000
        ut1, ut2 = 0.000, 0.000
    elif n == 2:
        ux1, ux2 = 0.500, 0.000
        ut1, ut2 = 0.000, 0.300
    elif n == 3:
        ux1, ux2 = 0.500, 0.000
        ut1, ut2 = 0.000, 0.500
    elif n == 4:
        ux1, ux2 = 0.500, 0.000
        ut1, ut2 = 0.000, 0.700
    elif n == 5:
        ux1, ux2 = 0.500, 0.000
        ut1, ut2 = 0.000, 0.900
    elif n == 6:
        ux1, ux2 = 0.500, 0.000
        ut1, ut2 = 0.000, 0.990
    elif n == 7:
        ux1, ux2 = 0.500, 0.000
        ut1, ut2 = 0.000, 0.999


    elif n == 8:
        ux1, ux2 = 0.000, 0.500
        ut1, ut2 = 0.000, 0.000
    elif n == 9:
        ux1, ux2 = 0.000, 0.500
        ut1, ut2 = 0.300, 0.000
    elif n == 10:
        ux1, ux2 = 0.000, 0.500
        ut1, ut2 = 0.500, 0.000
    elif n == 11:
        ux1, ux2 = 0.000, 0.500
        ut1, ut2 = 0.700, 0.000
    elif n == 12:
        ux1, ux2 = 0.000, 0.500
        ut1, ut2 = 0.900, 0.000
    elif n == 13:
        ux1, ux2 = 0.000, 0.500
        ut1, ut2 = 0.990, 0.000
    elif n == 14:
        ux1, ux2 = 0.000, 0.500
        ut1, ut2 = 0.999, 0.000

    else:
        ux1, ux2 = 0.000, 0.000
        ut1, ut2 = 0.000, 0.000

    if invert:
        rho1, ux1, ut1, p1, rho2, ux2, ut2, p2 = rho2, -ux2, ut2, p2, rho1, -ux1, ut1, p1

    solver = eval(rho1, ux1, ut1, p1, rho2, ux2, ut2, p2, gamma, T)
    print("%.3f\t%.3f\t%.3f\t%.3f" %(solver.solution.p3, solver.solution.ux3, solver.solution.rho3, solver.solution.rho4))


def marti_problem(n, T=0.4, invert=False, eval=simple_eval):
    print('marti ', n)
    gamma = 5. / 3

    rho1 = 1
    p1 = 1000
    p2 = 0.01
    rho2 = 1
    ux1 = 0
    ux2 = 0

    if n == 1:
        ut1, ut2 = 0.000, 0.000
    elif n == 2:
        ut1, ut2 = 0.000, 0.900
    elif n == 3:
        ut1, ut2 = 0.000, 0.990

    elif n == 4:
        ut1, ut2 = 0.900, 0.000
    elif n == 5:
        ut1, ut2 = 0.900, 0.900
    elif n == 6:
        ut1, ut2 = 0.900, 0.990

    elif n == 7:
        ut1, ut2 = 0.990, 0.000
    elif n == 8:
        ut1, ut2 = 0.990, 0.900
    elif n == 9:
        ut1, ut2 = 0.990, 0.990


    else:
        ux1, ux2 = 0.000, 0.000
        ut1, ut2 = 0.000, 0.000

    if invert:
        rho1, ux1, ut1, p1, rho2, ux2, ut2, p2 = rho2, -ux2, ut2, p2, rho1, -ux1, ut1, p1

    solver = eval(rho1, ux1, ut1, p1, rho2, ux2, ut2, p2, gamma, T)
    print("%.2e\t%.2e\t%.2e\t%.3f\t%.3f\t%.3f\t%.3f" %
          (solver.solution.rho3, solver.solution.rho4, solver.solution.p3, solver.solution.ux3,
          solver.solution.region_interfaces[3], solver.solution.region_interfaces[0], solver.solution.region_interfaces[1]))




if __name__ == '__main__':
    eval=simple_eval
    eval=complete_eval
    for i in range(5, 6):
        print('-------------------------------')
        #rezzolla_problem(i, eval=eval)
        marti_problem(i, eval=eval)

        print('-------------------------------\n')

        plt.show()

