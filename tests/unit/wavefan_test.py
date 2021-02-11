import srrp
import unittest
import numpy as np


class WavefanTest(unittest.TestCase):
    def test_wavefan_state_access(self):
        states = [srrp.State(0, 0, 0, 0), srrp.State(0.2, 0.2, 0.2, 0.2), srrp.State(0.5, 0.5, 0.5, 0.5),
                  srrp.State(0.7, 0.3, 0.3, 0.6)]
        waves = [srrp.ContactDiscontinuity(-0.3, 1), srrp.ContactDiscontinuity(0.1, 1),
                 srrp.ContactDiscontinuity(0.5, 1)]
        wavefan = srrp.Wavefan(states, waves)

        xis = np.array([-0.6, 0, 0.2, 0.7])

        self.assertEqual(wavefan.getRegionIndex(xis[0]), 0)
        self.assertEqual(wavefan.getRegionIndex(xis[1]), 1)
        self.assertEqual(wavefan.getRegionIndex(xis[2]), 2)
        self.assertEqual(wavefan.getRegionIndex(xis[3]), 3)

        for i in range(4):
            self.assertEqual(wavefan.states[i], states[i])

        for i in range(4):
            self.assertEqual(wavefan.regionStates[i](0), states[i])

        states2 = wavefan.getState(xis)

        def checkState(idx):
            self.assertEqual(states[idx].rho, states2.rho[idx])
            self.assertEqual(states[idx].vx, states2.vx[idx])
            self.assertEqual(states[idx].vt, states2.vt[idx])
            self.assertEqual(states[idx].pressure, states2.pressure[idx])

        checkState(0)
        checkState(1)
        checkState(2)
        checkState(3)


if __name__ == '__main__':
    unittest.main()
