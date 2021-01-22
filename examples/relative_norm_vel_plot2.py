from srrp import SRHDRiemannSolver
import matplotlib.pyplot as plt
import numpy as np


def get_dus(rhoL, uxL, utL, pL, rhoR, uxR, utR, pR, gamma, ps=[]):
    solver = SRHDRiemannSolver(rhoL, uxL, utL, pL, rhoR, uxR, utR, pR, gamma, 0, verbose=True)
    dus = []
    for p in ps:
        if p < pL:
            dus.append(solver.get_du_RR(p))
        elif p <= pR:
            dus.append(solver.get_du_RS(p))
        else:
            dus.append(solver.get_du_SS(p))

    return dus, solver.get_du_limit_RS(), solver.get_du_limit_SS()


lineArtist = lambda c='k': plt.Line2D((0,1),(0,0), color=c)

if __name__ == '__main__':
    rhoL = 1
    pL = 1
    pR = 0.1
    rhoR = 0.125

    #uxL = 0.
    #uxR = 0.3

    utR = 0
    gamma = 5/3



    for utL, uxL in [[0.00, 0.00], [0.90,0.0], [0.0, 0.9]]:
        ps = []
        dus = []
        for uxR in np.linspace(-0.99, 1, 100, endpoint=False):
            solver = SRHDRiemannSolver(rhoL, uxL, utL, pL, rhoR, uxR, utR, pR, gamma, 0, verbose=False)
            du = (uxL - uxR) / (1 - uxL * uxR)
            p = solver.solution.p3

            ps.append(p)
            dus.append(du)


        plt.plot(ps, dus, zorder=1, label="$v^t_1=$%.2f | $v^x_1=$%.2f"%(utL, uxL))
        #plt.scatter([pL, pR], [duRS, duRR], zorder=2, s=40, c='k')

    plt.ylim(-1.1, 1.1)
    plt.xlim(-0.1, 2.1)

    plt.plot([pL, pL], [-2,2], 'k--')
    plt.plot([pR, pR], [-2,2], 'k--')

    plt.grid(True)
    plt.legend(loc="lower right")

    plt.savefig('relative_velocities2.pdf')

    plt.show()


