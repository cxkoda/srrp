import numpy as np


class IdealEquationOfState:
    def __init__(self, gamma):
        self.gamma = gamma
        self.sigma = gamma / (gamma - 1)

    def computeSpeedOfSound(self, pressure, rho):
        h = self.computeEnthalpy(pressure, rho)
        return np.sqrt(self.gamma * pressure / (h * rho))

    def computeEnthalpy(self, pressure, rho):
        return 1 + self.sigma * pressure / rho

    def computeRho(self, pressure, h):
        return self.sigma / (h - 1) * pressure


class PolytropicEquationOfState:
    @classmethod
    def fromState(cls, state, gamma):
        K = PolytropicEquationOfState.computeIsentropicConstant(state, gamma)
        return cls(K, gamma)

    @staticmethod
    def computeIsentropicConstant(state, gamma):
        return state.pressure / np.power(state.rho, gamma)

    def __init__(self, isentropicConstant, gamma):
        self.isentropicConstant = isentropicConstant
        self.gamma = gamma
        self.sigma = gamma / (gamma - 1)

    def computeRho(self, pressure):
        return np.power(pressure / self.isentropicConstant, 1. / self.gamma)

    def computeSpeedOfSound(self, pressure, rho=None, h=None):
        if rho is None:
            rho = self.computeRho(pressure)
        if h is None:
            h = self.computeEnthalpy(pressure, rho)
        return np.sqrt(self.gamma * pressure / (h * rho))

    def computeEnthalpy(self, pressure, rho=None):
        if rho is None:
            rho = self.computeRho(pressure)
        return 1 + self.sigma * pressure / rho
