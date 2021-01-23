import numpy as np


class Shock:
    def computeJ(self, stateA, stateB):
        return np.sqrt(
            (stateA.pressure - stateB.pressure) / (self.eos.h(stateA) / stateA.rho - self.eos.h(stateB) / stateB.rhoB))

    def computeShockSpeed(self, stateA, J, sign):
        D = stateA.rho * stateA.lorentz()

        return (D ** 2 * stateA.vx + sign * J * np.sqrt(J ** 2 + D ** 2 * (1 - stateA.vx ** 2))) / (D ** 2 + J ** 2)

    def __init__(self, stateA, stateB, eos, sign):
        self.eos = eos
        j = self.computeJ(stateA, stateB)
        self.speed = self.computeShockSpeed(stateA, j, sign)
