import matplotlib.pyplot as plt
import numpy as np
import srrp

'''
Reproduction of (parts of) Figure 5 from
Rezzolla, Zanott, Pons (2003), J.Fluid Mech. 479, 199
DOI: 10.1017/S0022112002003506
'''

if __name__ == '__main__':
    gamma = 5 / 3
    stateL = srrp.State(rho=1, pressure=1, vx=0.0, vt=0.999)
    stateR = srrp.State(rho=0.125, pressure=0.1, vx=0.5, vt=0)
    xs = np.linspace(0, 1, 500)
    t = 0.4
    x0 = 0.5

    solver = srrp.Solver()
    solution = solver.solve(stateL, stateR, gamma)

    for wave in solution.waves:
        print(wave)

    for idx, state in enumerate(solution.states):
        print(f'{idx}: {state}')

    xis = (xs - x0) / t
    states = solution.getState(xis)
    plt.plot(xs, states.rho, label="rho")
    plt.plot(xs, states.pressure, label="pressure")
    plt.plot(xs, states.vx, label="vx")
    plt.plot(xs, states.vt, label="vt")

    plt.legend(loc='upper right')
    plt.grid(True, ls='--')
    plt.xlabel('x')

    plt.show()
