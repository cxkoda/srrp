# Special Relativistic Riemann Problem Solver
[![PyPI version shields.io](https://img.shields.io/pypi/v/srrp.svg)](https://pypi.org/project/srrp)
[![PyPI license](https://img.shields.io/pypi/l/srrp.svg)](https://pypi.org/project/srrp)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/srrp.svg)](https://pypi.org/project/srrp)
[![Test](https://github.com/cxkoda/srrp/workflows/Test/badge.svg)](https://github.com/cxkoda/srrp/actions)


One-dimensional Riemann problems have been proven to be useful tools to benchmark numerical hydrodynamic codes since they can be solved to arbitrary precision.
This library implements such a solver for special relativistic hydrodynamics based on the scheme presented by Rezzolla et al. (2003) [[doi: 10.1017/S0022112002003506](https://doi.org/10.1017/S0022112002003506)].

### Installation

This package can be obtained via pypi
```
pip install srrp
```


### Example Usage
The solution of a Riemann Problem producing a rarefaction-shock pattern could look like the following code (see also `examples/rarefaction-shock.py`)

```python
import matplotlib.pyplot as plt
import numpy as np
import srrp

'''
Reproduction of (parts of) Figure 3 from
Rezzolla, Zanott, Pons (2003), J.Fluid Mech. 479, 199
DOI: 10.1017/S0022112002003506
'''

if __name__ == '__main__':
    gamma = 5/3
    stateL = srrp.State(rho=1, pressure=1, vx=0.5, vt=0)
    stateR = srrp.State(rho=0.125, pressure=0.1, vx=0, vt=0)
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
```

For more examples visit the `examples` folder.
