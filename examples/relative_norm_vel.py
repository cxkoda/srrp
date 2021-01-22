from srrp import SRHDRiemannSolver
import matplotlib.pyplot as plt
import numpy as np


def plotDU(rhoL, uxL, utL, pL, rhoR, uxR, utR, pR, gamma, ps=[]):
    solver = SRHDRiemannSolver(rhoL, uxL, utL, pL, rhoR, uxR, utR, pR, gamma, 0, verbose=True)



    for lab, f in zip(['RR', 'RS', 'SS'], [solver.get_du_RS, solver.get_du_RR, solver.get_du_SS]):
        plt.plot(ps, [f(p) for p in ps], label=lab)

    plt.legend(loc=4)
    plt.ylim(-1, 1)

    vs = []
    for p in ps:
        if p < pR:
            vs.append(solver.get_du_RR(p))
        elif p < pL:
            vs.append(solver.get_du_RS(p))
        else:
            vs.append(solver.get_du_SS(p))

    plt.plot(ps, vs, 'k--')

    col_labels = ['L', 'R']
    row_labels = ['rho', 'ux', 'ut', 'p']
    table_vals = [[rhoL, rhoR], [uxL, uxR], [utL, utR], [pL, pR]]

    plt.table(cellText=table_vals,
              colWidths=[0.1] * 3,
              rowLabels=row_labels,
              colLabels=col_labels,
              loc='lower center')


if __name__ == '__main__':

    ps = np.linspace(0.001, 2, 2000)

    for uxL in [0.0]:
        for utL, utR in [[0,0.9] , [0,0.99],[0.9, 0.999]]:
    #    for utL, utR in [[0,0]]:
            print('----------------')
            plt.figure()
            plotDU(1, uxL, utL, 1, 0.125, 0.0, utR, 0.1, 5./ 3, ps=ps)

    #plotDU(10, -0.5, 0, 20, 1, 0.6, 0, 10, 5. / 3)
    plt.show()


