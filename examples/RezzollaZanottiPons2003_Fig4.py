import matplotlib.pyplot as plt
import numpy as np
import srrp
from mpl_toolkits.mplot3d import Axes3D

'''
Reproduction of (parts of) Figure 4 from
Rezzolla, Zanott, Pons (2003), J.Fluid Mech. 479, 199
DOI: 10.1017/S0022112002003506
'''


def getSolution(stateL, stateR, xis, gamma=5 / 3):
    return srrp.Solver().solve(stateL, stateR, gamma).getState(xis)


if __name__ == '__main__':
    stateL = srrp.State(rho=1, pressure=1, vx=0.5, vt=0)
    stateR = srrp.State(rho=0.125, pressure=0.1, vx=0, vt=0)
    xs = np.linspace(0, 1, 200)
    xis = (xs - 0.5) / 0.4
    vtRs = np.linspace(0, 0.9, 100)

    pressures = np.array([getSolution(stateL, stateR, xis).pressure for stateR.vt in vtRs])

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    xx, yy = np.meshgrid(xs, vtRs)
    surf = ax.plot_wireframe(xx, yy, pressures, rstride=4, cstride=4)
    ax.set_xlabel("x")
    ax.set_ylabel("$v^t_R$")
    ax.set_zlabel("p")

    plt.show()
