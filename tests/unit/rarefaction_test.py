import srrp
import unittest
import numpy as np
from srrp.EquationOfState import IdealEquationOfState


class RarefactionTest(unittest.TestCase):
    def test_state_inside_rarefaction(self):
        stateA = srrp.State(rho=1, vx=0, vt=0, pressure=10)
        eos = IdealEquationOfState(5 / 3)
        rarefaction = srrp.Rarefaction.fromStateAheadAndSpeedPressureBehind(stateA, 0, 1, eos, sign=+1)
        xis = np.linspace(rarefaction.speedTail, rarefaction.speedHead, 10)

        for xi in xis:
            state = rarefaction.computeRarefactionState(xi)
            self.assertLessEqual(state.rho, stateA.rho)
            self.assertLessEqual(state.pressure, stateA.pressure)


if __name__ == '__main__':
    unittest.main()
