import parameterized
import unittest
import srrp


@parameterized.parameterized_class(
    ('tangentialSpeedL', 'tangentialSpeedR', 'rhoLPrime', 'rhoRPrime', 'contactDiscontinuityPressure',
     'contactDiscontinuitySpeed', 'shockSpeed', 'rarefactionHeadSpeed', 'rarefactionTailSpeed'), [
        # (0.00, 0.00, 9.16e-2, 1.04e+1, 1.86e+1, 0.960, 0.987, -0.816, +0.668),
        # (0.00, 0.90, 1.51e-1, 1.46e+1, 4.28e+1, 0.913, 0.973, -0.816, +0.379),
        # (0.00, 0.99, 2.89e-1, 2.36e+1, 1.27e+2, 0.767, 0.927, -0.816, -0.132),
        (0.90, 0.00, 5.83e-3, 3.44e+0, 1.89e-1, 0.328, 0.452, -0.525, +0.308),
        # (0.90, 0.90, 1.49e-2, 4.46e+0, 9.04e-1, 0.319, 0.445, -0.525, +0.282),
        # (0.90, 0.99, 5.72e-2, 7.83e+0, 8.48e+0, 0.292, 0.484, -0.525, +0.197),
        # (0.99, 0.00, 1.99e-3, 1.91e+0, 3.16e-2, 0.099, 0.208, -0.196, +0.096),
        # (0.99, 0.90, 3.80e-3, 2.90e+0, 9.27e-2, 0.098, 0.153, -0.196, +0.094),
        # (0.99, 0.99, 1.29e-2, 4.29e+0, 7.06e-1, 0.095, 0.140, -0.196, +0.085),
    ])
class PonsTest(unittest.TestCase):
    """
    Testing against the tabulated solutions for Riemann problems in Pons et al. (2000) - Table 1
    doi: 10.1017/S0022112000001439
    The third row has probably a typo rhoRPrime = 43.6 -> 23.6
    """
    pressureL = 1e3
    rhoL = 1
    normalSpeedL = 0
    pressureR = 1e-2
    rhoR = 1
    normalSpeedR = 0
    gamma = 5 / 3

    def getSolution(self):
        if not hasattr(self, 'solution'):
            stateL = srrp.State(rho=self.rhoL, vx=self.normalSpeedL, vt=self.tangentialSpeedL,
                                pressure=self.pressureL)
            stateR = srrp.State(rho=self.rhoR, vx=self.normalSpeedR, vt=self.tangentialSpeedR,
                                pressure=self.pressureR)
            self.solution = srrp.Solver().solve(stateL, stateR, self.gamma)

        return self.solution

    def assertMyAlmostEqual(self, expected, actual):
        self.assertAlmostEqual(expected, actual, delta=abs(expected * 1e-2))

    def test_primedLeft(self):
        solution = self.getSolution()
        self.assertMyAlmostEqual(self.rhoLPrime, solution.states[1].rho)
        self.assertMyAlmostEqual(self.contactDiscontinuityPressure, solution.states[1].pressure)
        self.assertMyAlmostEqual(self.contactDiscontinuitySpeed, solution.states[1].vx)

    def test_primedRight(self):
        solution = self.getSolution()
        self.assertMyAlmostEqual(self.rhoRPrime, solution.states[2].rho)
        self.assertMyAlmostEqual(self.contactDiscontinuityPressure, solution.states[2].pressure)
        self.assertMyAlmostEqual(self.contactDiscontinuitySpeed, solution.states[2].vx)

    def test_rarefaction(self):
        solution = self.getSolution()
        self.assertTrue(isinstance(solution.waves[0], srrp.Rarefaction))
        self.assertMyAlmostEqual(self.rarefactionHeadSpeed, solution.waves[0].speedHead)
        self.assertMyAlmostEqual(self.rarefactionTailSpeed, solution.waves[0].speedTail)

    def test_shock(self):
        solution = self.getSolution()
        self.assertTrue(isinstance(solution.waves[2], srrp.Shock))
        self.assertMyAlmostEqual(self.shockSpeed, solution.waves[2].speed)

    def test_contactDiscontinuity(self):
        solution = self.getSolution()
        self.assertTrue(isinstance(solution.waves[1], srrp.ContactDiscontinuity))
        self.assertMyAlmostEqual(self.contactDiscontinuitySpeed, solution.waves[1].speed)


if __name__ == '__main__':
    unittest.main()
