import numpy as np


class State:
    def __init__(self, rho, vx, vt, pressure):
        self.rho = rho
        self.vx = vx
        self.vt = vt
        self.pressure = pressure

    def speed(self):
        return np.sqrt(self.vx ** 2 + self.vt ** 2)