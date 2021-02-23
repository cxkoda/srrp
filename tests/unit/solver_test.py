import srrp
import unittest
import numpy as np


class SolverTest(unittest.TestCase):
    def test_SameLeftRightState(self):
        stateL = srrp.State(rho=1, vx=0.5, vt=0, pressure=1)
        stateR = srrp.State(rho=0.125, vx=0, vt=0, pressure=0.1)
        gamma = 5 / 3
        solution = srrp.Solver().solve(stateL, stateR, gamma)

        self.assertEqual(solution.states[0], stateL)
        self.assertEqual(solution.states[-1], stateR)

    def test_reversed(self):
        stateL = srrp.State(rho=1, vx=0.5, vt=0, pressure=1)
        stateR = srrp.State(rho=0.125, vx=0, vt=0, pressure=0.1)
        gamma = 5 / 3
        solution1 = srrp.Solver().solve(stateL, stateR, gamma)

        stateL.vx *= -1
        stateR.vx *= -1
        solution2 = srrp.Solver().solve(stateR, stateL, gamma)

        for state1, state2 in zip(solution1.states, reversed(solution2.states)):
            self.assertEqual(state1.rho, state2.rho)
            self.assertEqual(state1.pressure, state2.pressure)
            self.assertEqual(state1.vt, state2.vt)
            self.assertEqual(state1.vx, -state2.vx)

        self.assertEqual(solution1.waves[0].speedHead, -solution2.waves[2].speedHead)
        self.assertEqual(solution1.waves[0].speedTail, -solution2.waves[2].speedTail)
        self.assertEqual(solution1.waves[1].speed, -solution2.waves[1].speed)
        self.assertEqual(solution1.waves[1].speed, -solution2.waves[1].speed)
        self.assertEqual(solution1.waves[2].speed, -solution2.waves[0].speed)


if __name__ == '__main__':
    unittest.main()
