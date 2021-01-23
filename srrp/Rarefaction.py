import numpy as np


class Shock:
    def rho_isentrope(self, p, K):
        return (p / K) ** (1. / self.gamma)

    def get_isentrope_K(self, rho, p):
        return p / (rho ** self.gamma)

    def cs_isentrope(self, p, K):
        ret = np.sqrt((self.gamma - 1) /
                       (1. / (K * self.sigma) *
                        (p / K) ** ((1 - self.gamma) / self.gamma) + 1))
        return ret

    def get_A(self, hA, WA, utA):
        return hA*WA*utA

    def ux(self, xi, p, K, A, sign=+1):
        cs = self.cs_isentrope(p, K)
        h = self.h(self.rho_isentrope(p, K), p)
        a = cs * h
        b = sign * np.sqrt(A**2 * (1 - cs**2) + h**2)
        return (a - b * xi) / (a * xi - b)

    def uxb_raref(self, uxA, pA, pB, KA, AA, sign=+1):
        A_sqr = AA*AA

        def integrand(p):
            rho = self.rho_isentrope(p, KA)
            h = self.h(rho, p)
            c = self.cs(rho, p)
            h_sqr = h*h
            c_sqr = c*c
            return np.sqrt(h_sqr + A_sqr * (1 - c_sqr)) / ((h_sqr + A_sqr) * rho * c)

        B = 0.5 * np.log((1 + uxA)/(1 - uxA)) + sign * integrate.quad(integrand, pA, pB)[0]
        return np.tanh(B)

    def xi_interface(self, ux, ut, c, sign=+1):
        u_sqr = ut*ut + ux*ux
        return (ux * (1 - c ** 2) + sign * c * np.sqrt(
            (1 - u_sqr) * (1 - u_sqr * c ** 2 - ux ** 2 * (1 - c ** 2)))) / (1 - u_sqr * c ** 2)

    def __init__(self, stateA, stateB, eos, sign):
        self.eos = eos
        j = self.computeJ(stateA, stateB)
        self.speed = self.computeShockSpeed(stateA, j, sign)


        self.K1 = self.get_isentrope_K(self.rho1, self.p1)
        self.c1 = self.cs(self.rho1, self.p1)

        self.rho3 = self.rho_isentrope(self.p3, self.K1)
        self.h3 = self.h(self.rho3, self.p3)

        self.h4 = self.get_h_taub(self.rho6, self.p6, self.h6, self.p4)
        self.rho4 = self.rho_h(self.p4, self.h4)

        self.ut3 = self.get_utb(self.A1, self.ux3, self.h3)
        self.ut4 = self.get_utb(self.A6, self.ux4, self.h4)

        self.c3 = self.cs(self.rho3, self.p3)
        vr2H = self.xi_interface(self.ux1, self.ut1, self.c1, sign=-1)
        vr2T = self.xi_interface(self.ux3, self.ut3, self.c3, sign=-1)

