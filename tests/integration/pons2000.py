import parameterized
import unittest
import srrp


@parameterized.parameterized_class(
    ('tangentialSpeedL', 'tangentialSpeedR', 'rhoLPrime', 'rhoRPrime', 'contactDiscontinuityPressure',
     'contactDiscontinuitySpeed', 'shockSpeed', 'rarefactionHeadSpeed', 'rarefactionTailSpeed'), [
        (0.00, 0.00, 9.16e-2, 1.04e+1, 1.86e+1, 0.960, 0.987, -0.816, +0.668),
        (0.00, 0.90, 1.51e-1, 1.46e+1, 4.28e+1, 0.913, 0.973, -0.816, +0.379),
        (0.00, 0.99, 2.89e-1, 2.36e+1, 1.27e+2, 0.767, 0.927, -0.816, -0.132),
        (0.90, 0.00, 5.83e-3, 3.44e+0, 1.89e-1, 0.328, 0.452, -0.525, +0.308),
        (0.90, 0.90, 1.49e-2, 4.46e+0, 9.04e-1, 0.319, 0.445, -0.525, +0.282),
        (0.90, 0.99, 5.72e-2, 7.83e+0, 8.48e+0, 0.292, 0.484, -0.525, +0.197),
        (0.99, 0.00, 1.99e-3, 1.91e+0, 3.16e-2, 0.099, 0.208, -0.196, +0.096),
        (0.99, 0.90, 3.80e-3, 2.90e+0, 9.27e-2, 0.098, 0.153, -0.196, +0.094),
        (0.99, 0.99, 1.29e-2, 4.29e+0, 7.06e-1, 0.095, 0.140, -0.196, +0.085),
    ])
class PonsTest(unittest.TestCase):
    """
    Testing against the tabulated solutions for Riemann problems in Pons et al. (2000) - Table 1
    doi: 10.1017/S0022112000001439
    The third row has probably a typo rhoRPrime = 43.6 -> 23.6
    """
    pressureL = 1e3
    rhoL = 1
    vxL = 0
    pressureR = 1e-2
    rhoR = 1
    vxR = 0
    gamma = 5 / 3

    def getSolution(self):
        if not hasattr(self, 'solver'):
            self.solver = srrp.SRHDRiemannSolver(self.rhoL, self.vxL, self.tangentialSpeedL, self.pressureL, self.rhoR,
                                                 self.vxR,
                                                 self.tangentialSpeedR,
                                                 self.pressureR, self.gamma)
        return self.solver.solution

    def assertMyAlmostEqual(self, expected, actual):
        self.assertAlmostEqual(expected, actual, delta=abs(expected * 1e-2))

    def test_primedLeft(self):
        solution = self.getSolution()
        self.assertMyAlmostEqual(self.rhoLPrime, solution.rho3)
        self.assertMyAlmostEqual(self.contactDiscontinuityPressure, solution.p3)
        self.assertMyAlmostEqual(self.contactDiscontinuitySpeed, solution.ux3)

    def test_primedRight(self):
        solution = self.getSolution()
        self.assertMyAlmostEqual(self.rhoRPrime, solution.rho4)
        self.assertMyAlmostEqual(self.contactDiscontinuityPressure, solution.p4)
        self.assertMyAlmostEqual(self.contactDiscontinuitySpeed, solution.ux4)

    def test_rarefaction(self):
        solution = self.getSolution()
        self.assertMyAlmostEqual(self.rarefactionHeadSpeed, solution.region_interfaces[0])
        self.assertMyAlmostEqual(self.rarefactionTailSpeed, solution.region_interfaces[1])

    def test_shock(self):
        solution = self.getSolution()
        self.assertMyAlmostEqual(self.shockSpeed, solution.region_interfaces[3])

    def test_contactDiscontinuity(self):
        solution = self.getSolution()
        self.assertMyAlmostEqual(self.contactDiscontinuitySpeed, solution.region_interfaces[2])


if __name__ == '__main__':
    unittest.main()
