import matplotlib.pyplot as plt
import numpy as np
import srrp

'''
Reproduction of (parts of) Figure 2 from
Rezzolla, Zanott, Pons (2003), J.Fluid Mech. 479, 199
DOI: 10.1017/S0022112002003506
'''

def get_dus(stateL, stateR, pressure, gamma=5/3):
    solver = srrp.Solver()
    solver.solve(stateL, stateR, gamma)
    dus = []

    for p in pressure:
        if p < stateR.pressure:
            dus.append(solver.get_du_RR(p))
        elif p < stateL.pressure:
            dus.append(solver.get_du_RS(p))
        else:
            dus.append(solver.get_du_SS(p))

    return dus, solver.get_du_limit_RS(), solver.get_du_limit_SS()

if __name__ == '__main__':
    stateL = srrp.State(rho=1, pressure=1, vx=0)
    stateR = srrp.State(rho=0.125, pressure=0.1, vx=.3)

    relativeNormalSpeed = srrp.Util.relativeSpeed(stateL.vx, stateR.vx)
    pressure = np.linspace(0, 2, 1000)

    for vtL, vtR in [[0.00, 0.00], [0.90,0.0], [0.99, 0.00],  [0.0, 0.90]]:
        stateL.vt = vtL
        stateR.vt = vtR
        dus, duRR, duRS = get_dus(stateL, stateR, pressure, gamma=5/3)

        plt.plot(pressure, dus, zorder=1, label="$v^t_1=$%.2f | $v^t_6=$%.2f"%(vtL, vtR))
        plt.scatter([stateL.pressure, stateR.pressure], [duRS, duRR], zorder=2, s=40, c='k')

    plt.ylim(-1.1, 1.1)
    plt.xlim(-0.1, 2.1)

    # Add helper lines
    plt.plot([-1, 3], [relativeNormalSpeed, relativeNormalSpeed], 'k--')
    plt.plot([stateL.pressure, stateL.pressure], [-2,2], 'k--')
    plt.plot([stateR.pressure, stateR.pressure], [-2,2], 'k--')

    plt.grid(True)
    plt.legend(loc="lower right")
    plt.xlabel("Contact Discontinuity Pressure")
    plt.ylabel("Left/Right - Relative Velocity")

    plt.show()


