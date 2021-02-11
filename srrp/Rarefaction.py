import numpy as np
from .State import State
from scipy import integrate, optimize
from .EquationOfState import PolytropicEquationOfState, IdealEquationOfState

from .Util import *


class Rarefaction:
    @classmethod
    def fromStateAheadAndSpeedPressureBehind(cls, stateA, speedB, pressureB, eos: IdealEquationOfState, sign):
        polytrope = PolytropicEquationOfState.fromState(stateA, eos.gamma)

        A = computeA(stateA, eos)
        h = polytrope.computeEnthalpy(pressureB)
        stateB = State()
        stateB.pressure = pressureB
        stateB.vx = speedB
        stateB.vt = computeVt(A, stateB.vx, h)
        stateB.rho = polytrope.computeRho(pressureB)
        return cls(stateA, stateB, eos, sign)

    @staticmethod
    def ux(xi, pressure, A, eos: PolytropicEquationOfState, sign):
        cs = eos.computeSpeedOfSound(pressure)
        h = eos.computeEnthalpy(pressure)
        a = cs * h
        b = sign * np.sqrt(A ** 2 * (1 - cs ** 2) + h ** 2)
        return (a - b * xi) / (a * xi - b)

    @staticmethod
    def computeVxb(stateA, pressureB, eos: PolytropicEquationOfState, sign, A=None):
        if A is None:
            A = computeA(stateA, eos)
        A_sqr = A ** 2

        def integrand(_pressure):
            rho = eos.computeRho(_pressure)
            h = eos.computeEnthalpy(_pressure, rho=rho)
            cs = eos.computeSpeedOfSound(_pressure, rho=rho, h=h)
            h_sqr = h ** 2
            cs_sqr = cs ** 2
            return np.sqrt(h_sqr + A_sqr * (1 - cs_sqr)) / ((h_sqr + A_sqr) * rho * cs)

        B = 0.5 * np.log((1 + stateA.vx) / (1 - stateA.vx)) + sign * \
            integrate.quad(integrand, stateA.pressure, pressureB)[0]
        return np.tanh(B)

    @staticmethod
    def xi_interface(state, eos: PolytropicEquationOfState, sign):
        speed_sqr = state.speed() ** 2
        cs = eos.computeSpeedOfSound(state.pressure)
        cs_sqr = cs ** 2
        return (state.vx * (1 - cs_sqr)
                + sign * cs * np.sqrt((1 - speed_sqr) * (1 - speed_sqr * cs_sqr - state.vx ** 2 * (1 - cs_sqr)))
                ) / (1 - speed_sqr * cs_sqr)

    def __init__(self, stateA, stateB, eos, sign):
        self.stateA = stateA
        self.stateB = stateB
        K = PolytropicEquationOfState.computeIsentropicConstant(stateA, eos.gamma)
        assert (K == PolytropicEquationOfState.computeIsentropicConstant(stateA, eos.gamma))

        self.eos = PolytropicEquationOfState(K, eos.gamma)
        self.sign = sign

        self.speedHead = self.xi_interface(stateA, self.eos, sign)
        self.speedTail = self.xi_interface(stateB, self.eos, sign)

        self.A = computeA(stateA, eos)

    def computeRarefactionState(self, xi):
        pmin = min(self.stateA.pressure, self.stateB.pressure)
        pmax = max(self.stateA.pressure, self.stateB.pressure)
        state = State()

        state.pressure = optimize.brentq(lambda pressure:
                                         self.ux(xi, pressure, self.A, self.eos, self.sign) -
                                         self.computeVxb(self.stateA, pressure, self.eos, self.sign, A=self.A),
                                         pmin, pmax)

        state.rho = self.eos.computeRho(state.pressure)
        state.vx = self.ux(xi, state.pressure, self.A, self.eos, self.sign)
        h = self.eos.computeEnthalpy(state.pressure, rho=state.rho)
        state.vt = computeVt(self.A, state.vx, h)
        return state

    def __str__(self):
        return f'Rarefaction: vHead={self.speedHead:.3f}, vTail={self.speedTail:.3f}'

    def __repr__(self):
        return str(self)
