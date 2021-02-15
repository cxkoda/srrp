import parameterized
import unittest
import srrp


@parameterized.parameterized_class(
    ('normalSpeedL', 'normalSpeedR', 'tangentialSpeedL', 'tangentialSpeedR', 'contactDiscontinuityPressure',
     'contactDiscontinuitySpeed', 'rhoLPrime', 'rhoRPrime', 'wavepattern'), [
        (0.5, 0.0, 0.0, 0.000, 0.597, 0.640, 0.734, 0.342, 'SR'),
        (0.5, 0.0, 0.0, 0.300, 0.621, 0.631, 0.751, 0.349, 'SR'),
        (0.5, 0.0, 0.0, 0.500, 0.673, 0.611, 0.788, 0.364, 'SR'),
        (0.5, 0.0, 0.0, 0.700, 0.787, 0.570, 0.866, 0.394, 'SR'),
        (0.5, 0.0, 0.0, 0.900, 1.150, 0.455, 1.088, 0.474, '2S'),
        (0.5, 0.0, 0.0, 0.990, 2.199, 0.212, 1.593, 0.647, '2S'),
        (0.5, 0.0, 0.0, 0.999, 3.011, 0.078, 1.905, 0.750, '2S'),
        (0.0, 0.5, 0.000, 0.0, 0.154, 0.620, 0.326, 0.162, 'SR'),
        (0.0, 0.5, 0.300, 0.0, 0.139, 0.594, 0.306, 0.152, 'SR'),
        (0.0, 0.5, 0.500, 0.0, 0.115, 0.542, 0.274, 0.136, 'SR'),
        (0.0, 0.5, 0.700, 0.0, 0.085, 0.450, 0.228, 0.113, '2R'),
        (0.0, 0.5, 0.900, 0.0, 0.051, 0.280, 0.168, 0.084, '2R'),
        (0.0, 0.5, 0.990, 0.0, 0.031, 0.095, 0.123, 0.061, '2R'),
        (0.0, 0.5, 0.999, 0.0, 0.026, 0.031, 0.110, 0.052, '2R'),
    ])
class RezzollaTest(unittest.TestCase):
    """
    Testing against the tabulated solutions for Riemann problems in Rezzolla et al. (2003) - Table 1
    doi: 10.1017/S0022112002003506
    """
    pressureL = 1
    rhoL = 1
    pressureR = 0.1
    rhoR = 0.125
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
        self.assertAlmostEqual(expected, actual, places=2)

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

    def test_wavepattern(self):
        solution = self.getSolution()
        if self.wavepattern == '2S':
            self.assertTrue(
                isinstance(solution.waves[0], srrp.Shock)
                and isinstance(solution.waves[-1], srrp.Shock)
            )
        elif self.wavepattern == '2R':
            self.assertTrue(
                isinstance(solution.waves[0], srrp.Rarefaction)
                and isinstance(solution.waves[-1], srrp.Rarefaction)
            )
        elif self.wavepattern == 'SR':
            self.assertTrue(
                (
                        isinstance(solution.waves[0], srrp.Shock)
                        and isinstance(solution.waves[-1], srrp.Rarefaction)
                ) or (
                        isinstance(solution.waves[0], srrp.Rarefaction)
                        and isinstance(solution.waves[-1], srrp.Shock)
                )
            )
        else:
            raise RuntimeError('No such wavepattern!')


if __name__ == '__main__':
    unittest.main()
