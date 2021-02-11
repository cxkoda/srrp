import scipy.optimize as opt
import scipy.integrate as integrate
import numpy as np

from .Shock import Shock
from .Rarefaction import Rarefaction
from .EquationOfState import IdealEquationOfState, PolytropicEquationOfState
from .Util import *

from .ContactDiscontinuity import ContactDiscontinuity
from .State import State
from .Wavefan import Wavefan


def getWavefan(state1, state6, cdSpeed, cdPressure, eos: IdealEquationOfState, waveLType, waveRtype, reversed=False):
    waveL = waveLType.fromStateAheadAndSpeedPressureBehind(state1, cdSpeed, cdPressure, eos, sign=-1)
    state3 = waveL.stateB
    waveR = waveRtype.fromStateAheadAndSpeedPressureBehind(state6, cdSpeed, cdPressure, eos, sign=+1)
    state4 = waveR.stateB
    cd = ContactDiscontinuity(cdSpeed, cdPressure)
    states = [state1, state3, state4, state6]
    waves = [waveL, cd, waveR]

    return Wavefan(states, waves, reversed)


class Solver:
    '''
    Solve SRHD Riemann Problem with non-zero tangential velocity:
        - Determine wave pattern (Shock-Shock, Rarefaction-Shock, etc.)
        - Compute p_star (pressure of the contact discontinuity)

    see paper
        Rezzolla_Zanotti_2003_Fluid_Mech_479
        DOI: 10.1017/S0022112002003506
    '''

    def get_du_SS(self, p_star):
        ux3 = Shock.computeVxb(self.state1, p_star, self.eos, sign=-1)
        ux4 = Shock.computeVxb(self.state6, p_star, self.eos, sign=+1)
        v13 = relativeSpeed(self.state1.vx, ux3)
        v64 = relativeSpeed(self.state6.vx, ux4)

        return relativeSpeed(v13, v64)

    def get_du_RS(self, p_star):
        polytrope1 = PolytropicEquationOfState.fromState(self.state1, self.eos.gamma)
        ux3 = Rarefaction.computeVxb(self.state1, p_star, polytrope1, sign=-1)
        ux4 = Shock.computeVxb(self.state6, p_star, self.eos, sign=+1)
        v13 = relativeSpeed(self.state1.vx, ux3)
        v64 = relativeSpeed(self.state6.vx, ux4)
        return relativeSpeed(v13, v64)

    def get_du_RR(self, p_star):
        polytrope1 = PolytropicEquationOfState.fromState(self.state1, self.eos.gamma)
        polytrope6 = PolytropicEquationOfState.fromState(self.state6, self.eos.gamma)
        ux3 = Rarefaction.computeVxb(self.state1, p_star, polytrope1, sign=-1)
        ux4 = Rarefaction.computeVxb(self.state6, p_star, polytrope6, sign=+1)
        v13 = relativeSpeed(self.state1.vx, ux3)
        v64 = relativeSpeed(self.state6.vx, ux4)
        return relativeSpeed(v13, v64)

    def get_Vs_lim(self):
        '''
        see Appendix C in Paper
        '''
        h3p = Shock.computeTaubAdiabat(self.state6, self.state1.pressure, self.eos)
        h6 = self.eos.computeEnthalpy(self.state6.rho, self.state6.pressure)
        J23p_sqr = Shock.computeJ_sqr(self.state6.pressure, self.state1.pressure, h6, h3p, self.eos)
        J23p = np.sqrt(np.abs(J23p_sqr))
        return Shock.computeShockSpeed(self.state6, J23p, sign=+1)

    def get_du_limit_SS(self):
        Vs = self.get_Vs_lim()
        return (((self.state1.pressure - self.state6.pressure) * (1 - self.state6.vx * Vs)) /
                ((Vs - self.state6.vx) * (self.h6 * self.state6.rho * self.W6_sqr * (
                        1 - self.state6.vx ** 2) + self.state1.pressure - self.state6.pressure)))

    def get_du_limit_RS(self):
        # Rarefaction.computeVxb(self.state1, self.state6.pressure, self.eos)

        A1_sqr = computeA(self.state1, self.eos) ** 2
        polytrope1 = PolytropicEquationOfState.fromState(self.state1, self.eos.gamma)

        def integrand(_pressure):
            rho = polytrope1.computeRho(_pressure)
            h = polytrope1.computeEnthalpy(_pressure, rho=rho)
            cs = polytrope1.computeSpeedOfSound(_pressure, rho=rho, h=h)
            h_sqr = h ** 2
            cs_sqr = cs ** 2
            return np.sqrt(h_sqr + A1_sqr * (1 - cs_sqr)) / ((h_sqr + A1_sqr) * rho * cs)

        return np.tanh(
            integrate.quad(integrand, self.state1.pressure, self.state6.pressure)[0]
        )

    def get_du_limit_RR(self):
        A1_sqr = computeA(self.state1, self.eos) ** 2
        A6_sqr = computeA(self.state6, self.eos) ** 2
        polytrope1 = PolytropicEquationOfState.fromState(self.state1, self.eos.gamma)
        polytrope6 = PolytropicEquationOfState.fromState(self.state6, self.eos.gamma)

        def integrand1(_pressure):
            rho = polytrope1.computeRho(_pressure)
            h = polytrope1.computeEnthalpy(_pressure, rho=rho)
            cs = polytrope1.computeSpeedOfSound(_pressure, rho=rho, h=h)
            h_sqr = h ** 2
            cs_sqr = cs ** 2
            return np.sqrt(h_sqr + A1_sqr * (1 - cs_sqr)) / ((h_sqr + A1_sqr) * rho * cs)

        def integrand6(_pressure):
            rho = polytrope6.computeRho(_pressure)
            h = polytrope6.computeEnthalpy(_pressure, rho=rho)
            cs = polytrope6.computeSpeedOfSound(_pressure, rho=rho, h=h)
            h_sqr = h ** 2
            cs_sqr = cs ** 2
            return np.sqrt(h_sqr + A6_sqr * (1 - cs_sqr)) / ((h_sqr + A6_sqr) * rho * cs)

        v1c = np.tanh(integrate.quad(integrand1, self.state1.pressure, 0))[0]
        v6c = np.tanh(integrate.quad(integrand6, 0, self.state6.pressure))[0]
        return relativeSpeed(v1c, v6c)

    def solve(self, stateL, stateR, gamma):
        if stateL.pressure >= stateR.pressure:
            self.reversed = False
            self.state1 = stateL
            self.state6 = stateR
        else:
            self.reversed = True
            self.state1 = stateR
            self.state6 = stateL

        self.eos = IdealEquationOfState(gamma)

        self.K1 = PolytropicEquationOfState.computeIsentropicConstant(self.state1, gamma)
        self.K6 = PolytropicEquationOfState.computeIsentropicConstant(self.state1, gamma)
        self.W1 = self.state1.lorentz()
        self.W6 = self.state6.lorentz()
        self.W1_sqr = self.W1 ** 2
        self.W6_sqr = self.W6 ** 2

        self.c1 = self.eos.computeSpeedOfSound(self.state1.rho, self.state1.pressure)
        self.c6 = self.eos.computeSpeedOfSound(self.state6.rho, self.state6.pressure)
        self.h1 = self.eos.computeEnthalpy(self.state1.rho, self.state1.pressure)
        self.h6 = self.eos.computeEnthalpy(self.state6.rho, self.state6.pressure)

        return self.determine_wave_pattern()

    def determine_wave_pattern(self):
        du_0 = relativeSpeed(self.state1.vx, self.state6.vx)
        eps = 1e-15

        if du_0 <= self.get_du_limit_RR():
            self.solution_type = 'RR*'
            p_star = 0
            ux_star = 0
            solution = getWavefan(self.state1, self.state6, ux_star, p_star, Rarefaction, Rarefaction,
                                  reversed=self.reversed)

        elif du_0 <= self.get_du_limit_RS():
            self.solution_type = 'RR'
            p_min = (self.state6.pressure + eps) * eps
            p_max = self.state1.pressure
            p_star = opt.brentq(lambda p: self.get_du_RR(p) - du_0, p_min, p_max)
            polytrope = PolytropicEquationOfState.fromState(self.state6, self.eos.gamma)
            ux_star = Rarefaction.computeVxb(self.state6, p_star, polytrope, sign=+1)
            solution = getWavefan(self.state1, self.state6, ux_star, p_star, self.eos, Rarefaction, Rarefaction,
                                  reversed=self.reversed)

        elif du_0 <= self.get_du_limit_SS():
            self.solution_type = 'RS'
            p_min = self.state6.pressure + eps
            p_max = self.state1.pressure
            p_star = opt.brentq(lambda p: self.get_du_RS(p) - du_0, p_min, p_max)
            ux_star = Shock.computeVxb(self.state6, p_star, self.eos, sign=+1)
            solution = getWavefan(self.state1, self.state6, ux_star, p_star, self.eos, Rarefaction, Shock,
                                  reversed=self.reversed)

        else:
            self.solution_type = 'SS'
            p_star_guess = 1.1 * self.state1.pressure
            p_star = opt.root(lambda p: self.get_du_SS(p) - du_0, p_star_guess).x[0]
            ux_star = Shock.computeVxb(self.state6, p_star, self.eos, sign=+1)
            solution = getWavefan(self.state1, self.state6, ux_star, p_star, self.eos, Shock, Shock,
                                  reversed=self.reversed)

        return solution