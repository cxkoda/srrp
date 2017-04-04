from libSRHD import SRHDRiemannSolver
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm





def get_ps(rhoL, uxL, utL, pL, rhoR, uxR, utR, pR, gamma, T, xs, x0=0.5):
    print(rhoL, uxL, utL, pL, rhoR, uxR, utR, pR)
    solver = SRHDRiemannSolver(rhoL, uxL, utL, pL, rhoR, uxR, utR, pR, gamma, x0, verbose=False)
    qs = np.array([solver(x, T) for x in xs])
    return qs[:, 3]

def sod_problem2(N=100):
    gamma = 5./3
    rhoL = 1
    pL = 1
    uxL = 0.
    pR = 0.1
    rhoR = 0.125
    uxR = 0.5
    utL = None
    utR = 0
    variab = np.linspace(0.5, 0.999, N)
    return rhoL, uxL, utL, pL, rhoR, uxR, utR, pR, gamma, variab


def sod_problem(N=100):
    gamma = 5./3
    rhoL = 1
    pL = 1
    uxL = 0.5
    pR = 0.1
    rhoR = 0.125
    uxR = 0
    utL = 0
    utR = None
    variab = np.linspace(0, 0.9, N)
    return rhoL, uxL, utL, pL, rhoR, uxR, utR, pR, gamma, variab


if __name__ == '__main__':

    xs = np.linspace(0, 1, 100)
    T = 0.4

    rhoL, uxL, utL, pL, rhoR, uxR, utR, pR, gamma, variab = sod_problem2(50)
    args = [rhoL, uxL, utL, pL, rhoR, uxR, utR, pR, gamma]
    var_idx = [idx for idx, v in enumerate(args) if v is None][0]

    def set_arg(v):
        args[var_idx] = v
        return True

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ps = [get_ps(*args, T, xs, x0=0.5) for v in variab if set_arg(v)]
    xs, variab = np.meshgrid(xs, variab)
    surf = ax.plot_surface(xs, variab, ps, rstride=1, cstride=1, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)
    ax.set_xlabel('x')
    ax.set_ylabel('variable')
    ax.set_zlabel('p')

    col_labels = ['L', 'R']
    row_labels = ['rho', 'ux', 'ut', 'p']
    table_vals = [[rhoL, rhoR], [uxL, uxR], [utL, utR], [pL, pR]]

    plt.table(cellText=table_vals,
              colWidths=[0.1] * 3,
              rowLabels=row_labels,
              colLabels=col_labels,
              loc='lower left')

    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.show()


