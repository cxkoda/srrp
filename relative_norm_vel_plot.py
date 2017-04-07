from libSRHD import SRHDRiemannSolver
import matplotlib.pyplot as plt
import numpy as np


def get_dus(rhoL, uxL, utL, pL, rhoR, uxR, utR, pR, gamma, ps=[]):
    solver = SRHDRiemannSolver(rhoL, uxL, utL, pL, rhoR, uxR, utR, pR, gamma, 0, verbose=True)
    dus = []
    for p in ps:
        if p < pR:
            dus.append(solver.get_du_RR(p))
        elif p < pL:
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

    uxL = 0.
    uxR = 0.3

    du=(uxL - uxR)/(1-uxL*uxR)


    ps = np.linspace(0, 2, 1000)

    for utL, utR in [[0.00, 0.00], [0.90,0.0], [0.99, 0.00],  [0.0, 0.90]]:
        print('----------------')
        dus, duRR, duRS = get_dus(rhoL, uxL, utL, pL, rhoR, uxR, utR, pR, 5./ 3, ps=ps)

        plt.plot(ps, dus, zorder=1, label="$v^t_1=$%.2f | $v^t_6=$%.2f"%(utL, utR))
        plt.scatter([pL, pR], [duRS, duRR], zorder=2, s=40, c='k')

    plt.ylim(-1.1, 1.1)
    plt.xlim(-0.1, 2.1)

    plt.plot([-1, 3], [du, du], 'k--')
    plt.plot([pL, pL], [-2,2], 'k--')
    plt.plot([pR, pR], [-2,2], 'k--')

    plt.grid(True)
    plt.legend(loc="lower right")

    plt.show()


