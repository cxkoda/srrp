from libSRHD import SRHDRiemannSolver
import matplotlib.pyplot as plt
import numpy as np

def simple_eval(rhoL, uxL, utL, pL, rhoR, uxR, utR, pR, gamma, T, x0=0.5):
    solver = SRHDRiemannSolver(rhoL, uxL, utL, pL, rhoR, uxR, utR, pR, gamma, x0, verbose=True)
    #print(solver.solution.p3, solver.solution.u3, solver.solution.rho3, solver.solution.rho4)


def complete_eval(rhoL, uxL, utL, pL, rhoR, uxR, utR, pR, gamma, T, x0=0.5, dx=0.5, N=1000):
    solver = SRHDRiemannSolver(rhoL, uxL, utL, pL, rhoR, uxR, utR, pR, gamma, x0, verbose=True)
    xs = np.linspace(x0 - dx, x0 + dx, N)

    qs = np.array([solver(x, T) for x in xs])

    regIdxs = np.array([solver.getRegIdx(x, T) for x in xs])

    for i in range(1,7):
        sel = regIdxs == i
        for j in range(3):
            plt.figure(j)
            plt.plot(xs[sel], qs[sel][:, j])

    for j in range(3):
        plt.figure(j)
        plt.xlim(np.min(xs), np.max(xs))
        y_min, y_max = np.min(qs[:, j]), np.max(qs[:, j])
        dy = y_max - y_min
        plt.ylim(y_min - 0.1*dy, y_max + 0.1*dy)


def plotDU(rhoL, uxL, utL, pL, rhoR, uxR, utR, pR, gamma, ps=[]):
    solver = SRHDRiemannSolver(rhoL, uxL, utL, pL, rhoR, uxR, utR, pR, gamma, 0, verbose=True)



    for lab, f in zip(['RR', 'RS', 'SS'], [solver.get_du_RS, solver.get_du_RR, solver.get_du_SS]):
        plt.plot(ps, [f(p) for p in ps], label=lab)

    plt.legend(loc=4)
    plt.ylim(-1, 1)

    vs = []
    for p in ps:
        if p < pL:
            vs.append(solver.get_du_RR(p))
        elif p < pR:
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

    ps = np.linspace(0.001, 2, 1000)

    for uxL in [0.5]:
        for utL, utR in [[0,0.9] , [0,0.99],[0., 0.999]]:
    #    for utL, utR in [[0,0]]:
            print('----------------')
            plt.figure()
            plotDU(1, uxL, utL, 1, 0.125, 0.0, utR, 0.1, 5./ 3, ps=np.linspace(0.001, 2, 200))

    #plotDU(10, -0.5, 0, 20, 1, 0.6, 0, 10, 5. / 3)
    plt.show()


