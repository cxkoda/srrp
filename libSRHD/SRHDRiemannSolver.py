import scipy.optimize as opt
import scipy.integrate as integrate
import numpy as np

class SRHDRiemannSolver:
    '''
    Solve SRHD Riemann Problem with non-zero tangential velocity:
        - Determine wave pattern (Shock-Shock, Rarefaction-Shock, etc.)
        - Compute p_star (pressure of the contact discontinuity)

    see paper
        Rezzolla_Zanotti_2003_Fluid_Mech_479
        DOI: 10.1017/S0022112002003506
    '''
    def relative_u(self, uL, uR):
        return (uL - uR) / (1 - uL * uR)

    def W_sqr(self, u):
        return 1. / (1. - u*u)

    def W(self, u):
        return np.sqrt(self.W_sqr(u))

    def h(self, rho, p):
        return 1. + p/rho * self.sigma

    def cs(self, rho, p):
        return np.sqrt(self.gamma * p / (rho + self.sigma * p))


    def rho_isentrope(self, p, K):
        return np.power(p/K, 1./self.gamma)

    def get_isentrope_K(self, rho, p):
        return p / np.power(rho, self.gamma)

    def get_h_taub(self, rhoA, pA, hA, pB):
        c_2 = (1. + (pA - pB) / (pB * self.sigma))
        c_1 = - (pA - pB) / (pB * self.sigma)
        c_0 = hA * (pA - pB) / rhoA - hA ** 2
        if c_2 == 0:
            hB = - c_0/c_1
        else:
            hB = (-c_1 + np.sqrt(c_1*c_1 - 4*c_2*c_0)) / (2*c_2)

        return hB

    def get_utb(self, rhoA, uxA, utA, pA, WA, hA, pB, uxB, hB):
        return hA * WA * utA * np.sqrt((1 - uxB**2) / (hB**2 + (hA*WA*utA)**2))


    def get_J_sqr(self, rhoA, pA, hA, pB):
        hB = self.get_h_taub(rhoA, pA, hA, pB)
        return - self.sigma * (pA - pB) / (hA * (hA - 1.) / pA - hB * (hB - 1.) / pB)

    def get_Vs(self, rhoA, uxA, WA_sqr, J, J_sqr, sign=+1):
        return (rhoA**2 * WA_sqr * uxA + sign * J * np.sqrt(J_sqr + rhoA**2 * WA_sqr * (1 - uxA**2))) / (rhoA**2 * WA_sqr + J_sqr)


    def uxb_shock(self, rhoA, uxA, utA, pA, pB, WA, hA, sign=+1):
        J_sqr = self.get_J_sqr(rhoA, pA, hA, pB)
        J = np.sqrt(np.abs(J_sqr))
        Vs = self.get_Vs(rhoA, uxA, WA*WA, J, J_sqr, sign=sign)
        Ws = self.W(Vs)
        return (hA * WA * uxA + sign * Ws * (pB - pA)/J) / (hA * WA + (pB - pA)*(sign * Ws * uxA / J + 1/(rhoA * WA)))

    def uxb_raref(self, rhoA, uxA, utA, pA, pB, WA, hA, sign=+1):
        A = hA * WA * utA
        A_sqr = A*A

        K = self.get_isentrope_K(rhoA, pA)

        def integrand(p):
            rho = self.rho_isentrope(p, K)
            h = self.h(rho, p)
            c = self.cs(rho, p)
            h_sqr = h*h
            c_sqr = c*c
            return np.sqrt(h_sqr + A_sqr * (1 - c_sqr)) / ((h_sqr + A_sqr) * rho * c)

        B = 0.5 * np.log((1 + uxA)/(1 - uxA)) + sign * integrate.quad(integrand, pA, pB)[0]
        return np.tanh(B)


    def get_du_SS(self, p_star):
        ux3 = self.uxb_shock(self.rho1, self.ux1, self.ut1, self.p1, p_star, self.W1, self.h1, sign=-1)
        ux4 = self.uxb_shock(self.rho6, self.ux6, self.ut6, self.p6, p_star, self.W6, self.h6, sign=+1)
        v13 = self.relative_u(self.ux1, ux3)
        v64 = self.relative_u(self.ux6, ux4)

        return self.relative_u(v13, v64)

    def get_du_RS(self, p_star):
        ux3 = self.uxb_raref(self.rho1, self.ux1, self.ut1, self.p1, p_star, self.W1, self.h1, sign=-1)
        ux4 = self.uxb_shock(self.rho6, self.ux6, self.ut6, self.p6, p_star, self.W6, self.h6, sign=+1)
        v13 = self.relative_u(self.ux1, ux3)
        v64 = self.relative_u(self.ux6, ux4)
        return self.relative_u(v13, v64)

    def get_du_RR(self, p_star):
        ux3 = self.uxb_raref(self.rho1, self.ux1, self.ut1, self.p1, p_star, self.W1, self.h1, sign=-1)
        ux4 = self.uxb_raref(self.rho6, self.ux6, self.ut6, self.p6, p_star, self.W6, self.h6, sign=+1)
        v13 = self.relative_u(self.ux1, ux3)
        v64 = self.relative_u(self.ux6, ux4)
        return self.relative_u(v13, v64)


    def get_Vs_lim(self):
        '''
        see Appendix C in Paper
        '''
        h3p = self.get_h_taub(self.rho6, self.p6, self.h6, self.p1)
        J23p_sqr = self.sigma * (self.p6 - self.p1) / (h3p * (h3p-1)/self.p1 - self.h6*(self.h6-1)/self.p6)
        J23p = np.sqrt(np.abs(J23p_sqr))

        return self.get_Vs(self.rho6, self.ux6, self.W6_sqr, J23p, J23p_sqr, sign=+1)

    def get_du_limit_SS(self):

        Vs = self.get_Vs_lim()
        return (
            (
                (self.p1 - self.p6) * (1 - self.ux6 * Vs)
            ) / (
                (Vs - self.ux6) * (self.h6 * self.rho6 * self.W6_sqr * (1 - self.ux6**2) + self.p1 - self.p6)
            )
        )

    def get_du_limit_RS(self):
        A1_sqr = self.h1**2 * self.W1_sqr * self.ut1**2

        def integrand(p):
            rho = self.rho_isentrope(p, self.K1)
            h = self.h(rho, p)
            c = self.cs(rho, p)
            h_sqr = h*h
            c_sqr = c*c
            return np.sqrt(h_sqr + A1_sqr * (1 - c_sqr)) / ((h_sqr + A1_sqr) * rho * c)

        return np.tanh(
            integrate.quad(integrand, self.p1, self.p6)[0]
        )

    def get_du_limit_RR(self):
        A1_sqr = self.h1 ** 2 * self.W1_sqr * self.ut1 ** 2
        A6_sqr = self.h6 ** 2 * self.W6_sqr * self.ut6 ** 2

        def integrand1(p):
            rho = self.rho_isentrope(p, self.K1)
            h = self.h(rho, p)
            c = self.cs(rho, p)
            h_sqr = h*h
            c_sqr = c*c
            return np.sqrt(h_sqr + A1_sqr * (1 - c_sqr)) / ((h_sqr + A1_sqr) * rho * c)

        def integrand6(p):
            rho = self.rho_isentrope(p, self.K6)
            h = self.h(rho, p)
            c = self.cs(rho, p)
            h_sqr = h*h
            c_sqr = c*c
            return np.sqrt(h_sqr + A6_sqr * (1 - c_sqr)) / ((h_sqr + A6_sqr) * rho * c)

        v1c = np.tanh(integrate.quad(integrand1, self.p1, 0))
        v6c = np.tanh(integrate.quad(integrand6, 0, self.p6))
        return self.relative_u(v1c[0], v6c[0])

    def __init__(self, rhoL, uxL, utL, pL, rhoR, uxR, utR, pR, gamma, x0=0.5, verbose=False):
        self.solution_type = None
        uL = np.sqrt(uxL**2 + utL**2)
        uR = np.sqrt(uxR**2 + utR**2)


        if pL >= pR:
            self.reversed = False
            self.rho1 = rhoL
            self.rho6 = rhoR
            u1 = uL
            u6 = uR
            self.p1 = pL
            self.p6 = pR
            self.ux1 = uxL
            self.ut1 = utL
            self.ux6 = uxR
            self.ut6 = utR
        else:
            self.reversed = True
            self.rho6 = rhoL
            self.rho1 = rhoR
            u6 = uL
            u1 = uR
            self.p6 = pL
            self.p1 = pR
            self.ux1 = -uxR
            self.ut1 = utR
            self.ux6 = -uxL
            self.ut6 = utL

        self.verbose = verbose
        self.x0 = x0
        self.gamma = gamma
        self.sigma = self.gamma / (self.gamma - 1)

        self.K1 = self.get_isentrope_K(self.rho1, self.p1)
        self.K6 = self.get_isentrope_K(self.rho6, self.p6)
        self.W1_sqr = self.W_sqr(u1)
        self.W6_sqr = self.W_sqr(u6)
        self.W1 = np.sqrt(self.W1_sqr)
        self.W6 = np.sqrt(self.W6_sqr)

        self.c1 = self.cs(self.rho1, self.p1)
        self.c6 = self.cs(self.rho6, self.p6)
        self.h1 = self.h(self.rho1, self.p1)
        self.h6 = self.h(self.rho6, self.p6)

        self.determine_wave_pattern()

    def determine_wave_pattern(self):
        du_0 = self.relative_u(self.ux1, self.ux6)
        #print('du0: ', du_0, self.get_du_limit_RR(), self.get_du_limit_RS(), self.get_du_limit_SS())
        eps = 1e-15

        if du_0 <= self.get_du_limit_RR():
            self.solution_type = 'RR*'
            if self.verbose:
                print('Rarefaction-Rarefaction with Vacuum')
            p_star = 0
            ux_star = 0
            solution = RareRareSolution

        elif du_0 <= self.get_du_limit_RS():
            self.solution_type = 'RR'
            if self.verbose:
                print('Rarefaction-Rarefaction')
            p_min = (self.p6 + eps) * eps
            p_max = self.p1
            p_star = opt.brentq(lambda p: self.get_du_RR(p) - du_0, p_min, p_max)
            ux_star = self.uxb_raref(self.rho6, self.ux6, self.ut6, self.p6, p_star, self.W6, self.h6, sign=+1)
            solution = RareRareSolution

        elif du_0 <= self.get_du_limit_SS():
            self.solution_type = 'RS'
            if self.verbose:
                if self.reversed:
                    print('Shock-Rarefaction')
                else:
                    print('Rarefaction-Shock')
            p_min = self.p6 + eps
            p_max = self.p1
            p_star = opt.brentq(lambda p: self.get_du_RS(p) - du_0, p_min, p_max)
            ux_star = self.uxb_shock(self.rho6, self.ux6, self.ut6, self.p6, p_star, self.W6, self.h6, sign=+1)
            solution = RareShockSolution

        else:
            self.solution_type = 'SS'
            if self.verbose:
                print('Shock-Shock')
            p_star_guess = 1.1*self.p1
            p_star = opt.root(lambda p: self.get_du_SS(p) - du_0, p_star_guess).x[0]
            ux_star = self.uxb_shock(self.rho6, self.ux6, self.ut6, self.p6, p_star, self.W6, self.h6, sign=+1)
            solution = ShockShockSolution


        self.solution = solution(self.rho1, self.ux1, self.ut1, self.p1,
                                 ux_star, p_star,
                                 self.rho6, self.ux6, self.ut6, self.p6, self.gamma, self.x0)

    def __call__(self, x, t):
        if self.reversed:
            ret = self.solution(2 * self.x0 - x, t)
            ret[1] *= -1
            return ret
        else:
            return self.solution(x, t)

    def getRegIdx(self, x, t):
        if self.reversed:
            return self.solution.getRegIdx(2 * self.x0 - x, t)
        else:
            return self.solution.getRegIdx(x, t)



class SRHDSolution:
    '''
    Base-class for the construction of the actual waves with given p_star and ux_star
    '''


    # -------- Shock Relevant -----------

    def get_h_taub(self, rhoA, pA, hA, pB):
        c_2 = (1. + (pA - pB) / (pB * self.sigma))
        c_1 = - (pA - pB) / (pB * self.sigma)
        c_0 = hA * (pA - pB) / rhoA - hA ** 2
        if c_2 == 0:
            hB = - c_0 / c_1
        else:
            hB = (-c_1 + np.sqrt(c_1 * c_1 - 4 * c_2 * c_0)) / (2 * c_2)

        return hB

    def j(self, rhoA, pA, rhoB, pB):
        return np.sqrt((pB - pA) / (self.h(rhoA, pA)/rhoA - self.h(rhoB, pB)/rhoB))

    def get_Vs(self, rhoA, uxA, WA_sqr, J, J_sqr, sign=+1):
        return (rhoA**2 * WA_sqr * uxA + sign * J * np.sqrt(J_sqr + rhoA**2 * WA_sqr * (1 - uxA**2))) / (rhoA**2 * WA_sqr + J_sqr)

    def rho_h(self, p, h):
        return p * self.sigma / (h - 1)


    # ----------- Rarefaction relevant -----------

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

    # --------------- common --------------

    def W(self, u):
        return np.sqrt(1./(1 - u*u))

    def cs(self, rho, p):
        return np.sqrt(self.gamma * p / (rho + self.sigma * p))

    def h(self, rho, p):
        return 1. + p / rho * self.sigma

    def get_utb(self, AA, uxB, hB):
        return AA * np.sqrt((1 - uxB*uxB) / (hB*hB + AA*AA))

    def get_xi(self, x, t):
        return (x - self.x0) / t

    def region1(self, xi):
        return np.array([self.rho1, self.ux1, self.ut1, self.p1])

    def region2(self, xi):
        pmin = min(self.p1, self.p3)
        pmax = max(self.p1, self.p3)
        p2 = opt.brentq(lambda p:
                         self.ux(xi, p, self.K1, self.A1, sign=-1) -
                         self.uxb_raref(self.ux1, self.p1, p, self.K1, self.A1, sign=-1),
                      pmin, pmax)

        rho2 = self.rho_isentrope(p2, self.K1)
        ux2 = self.ux(xi, p2, self.K1, self.A1, sign=-1)
        h2 = self.h(rho2, p2)
        ut2 = self.get_utb(self.A1, ux2, h2)
        return np.array([rho2, ux2, ut2, p2])

    def region3(self, xi):
        return np.array([self.rho3, self.ux3, self.ut3, self.p3])

    def region4(self, xi):
        return np.array([self.rho4, self.ux4, self.ut4, self.p4])

    def region5(self, xi):
        pmin = min(self.p4, self.p6)
        pmax = max(self.p4, self.p6)
        p5 = opt.brentq(lambda p:
                         self.ux(xi, p, self.K6, self.A6, sign=+1) -
                         self.uxb_raref(self.ux6, self.p6, p, self.K6, self.A6, sign=+1),
                      pmin, pmax)


        rho5 = self.rho_isentrope(p5, self.K6)
        ux5 = self.ux(xi, p5, self.K6, self.A6, sign=+1)
        h5 = self.h(rho5, p5)
        ut5 = self.get_utb(self.A6, ux5, h5)
        return np.array([rho5, ux5, ut5, p5])

    def region6(self, xi):
        return np.array([self.rho6, self.ux6, self.ut6, self.p6])


    def __init__(self, rho1, ux1, ut1, p1, ux_star, p_star, rho6, ux6, ut6, p6, gamma, x0):
        self.rho1 = rho1
        self.rho6 = rho6
        self.ux1 = ux1
        self.ux6 = ux6
        self.ut1 = ut1
        self.ut6 = ut6
        self.p1 = p1
        self.p6 = p6
        self.u1 = np.sqrt(ux1 * ux1 + ut1 * ut1)
        self.u6 = np.sqrt(ux6 * ux6 + ut6 * ut6)

        self.p3 = p_star
        self.p4 = p_star
        self.ux3 = ux_star
        self.ux4 = ux_star

        self.x0 = x0
        self.gamma = gamma
        self.gamma_hat = np.sqrt(gamma - 1)
        self.sigma = self.gamma / (self.gamma - 1)

        self.W1 = self.W(self.u1)
        self.W6 = self.W(self.u6)

        self.h1 = self.h(self.rho1, self.p1)
        self.h6 = self.h(self.rho6, self.p6)

        self.A1 = self.get_A(self.h1, self.W1, self.ut1)
        self.A6 = self.get_A(self.h6, self.W6, self.ut6)

        self.region_interfaces = [None]*5
        self.regions = [self.region1, self.region2, self.region3, self.region4, self.region5]

    def __call__(self, x, t):
        xi = self.get_xi(x, t)

        for region, interface in zip(self.regions, self.region_interfaces):
            if interface is not None and xi <= interface:
                return region(xi)

        return self.region6(xi)

    def getRegIdx(self, x, t):
        xi = self.get_xi(x, t)

        for idx, interface in enumerate(self.region_interfaces):
            if interface is not None and xi <= interface:
                return idx + 1

        return 6

    def endInit(self):
        self.W3 = self.W(np.sqrt(self.ux3 ** 2 + self.ut3 ** 2))
        self.W4 = self.W(np.sqrt(self.ux4 ** 2 + self.ut4 ** 2))

class ShockShockSolution(SRHDSolution):
    def __init__(self, rho1, ux1, ut1, p1, ux_star, p_star, rho6, ux6, ut6, p6, gamma, x0):
        SRHDSolution.__init__(self, rho1, ux1, ut1, p1, ux_star, p_star, rho6, ux6, ut6, p6, gamma, x0)

        self.h3 = self.get_h_taub(self.rho1, self.p1, self.h1, self.p3)
        self.h4 = self.get_h_taub(self.rho6, self.p6, self.h6, self.p4)

        self.rho3 = self.rho_h(self.p3, self.h3)
        self.rho4 = self.rho_h(self.p4, self.h4)

        self.ut3 = self.get_utb(self.A1, self.ux3, self.h3)
        self.ut4 = self.get_utb(self.A6, self.ux4, self.h4)

        j2 = self.j(self.rho1, self.p1, self.rho3, self.p3)
        j5 = self.j(self.rho6, self.p6, self.rho4, self.p4)

        vs2 = self.get_Vs(self.rho1, self.ux1, self.W1**2, j2, j2**2, sign=-1)
        vs5 = self.get_Vs(self.rho6, self.ux6, self.W6**2, j5, j5**2, sign=+1)

        self.region_interfaces = [vs2, None, self.ux3, vs5, None]
        self.endInit()

class RareShockSolution(SRHDSolution):
    def __init__(self, rho1, ux1, ut1, p1, ux_star, p_star, rho6, ux6, ut6, p6, gamma, x0):
        SRHDSolution.__init__(self, rho1, ux1, ut1, p1, ux_star, p_star, rho6, ux6, ut6, p6, gamma, x0)

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

        j5 = self.j(self.rho6, self.p6, self.rho4, self.p4)
        vs5 = self.get_Vs(self.rho6, self.ux6, self.W6**2, j5, j5**2, sign=+1)

        self.region_interfaces = [vr2H, vr2T, self.ux3, vs5, None]
        self.endInit()

class RareRareSolution(SRHDSolution):
    def __init__(self, rho1, ux1, ut1, p1, ux_star, p_star, rho6, ux6, ut6, p6, gamma, x0):
        SRHDSolution.__init__(self, rho1, ux1, ut1, p1, ux_star, p_star, rho6, ux6, ut6, p6, gamma, x0)
        self.K1 = self.get_isentrope_K(self.rho1, self.p1)
        self.K6 = self.get_isentrope_K(self.rho6, self.p6)

        self.c1 = self.cs(self.rho1, self.p1)
        self.c6 = self.cs(self.rho6, self.p6)

        self.rho3 = self.rho_isentrope(self.p3, self.K1)
        self.rho4 = self.rho_isentrope(self.p4, self.K6)

        self.h3 = self.h(self.rho3, self.p3)
        self.h4 = self.h(self.rho4, self.p4)

        self.ut3 = self.get_utb(self.A1, self.ux3, self.h3)
        self.ut4 = self.get_utb(self.A6, self.ux4, self.h4)


        self.c3 = self.cs(self.rho3, self.p3)
        self.c4 = self.cs(self.rho4, self.p4)


        vr2H = self.xi_interface(self.ux1, self.ut1, self.c1, sign=-1)
        vr2T = self.xi_interface(self.ux3, self.ut3, self.c3, sign=-1)

        vr5H = self.xi_interface(self.ux6, self.ut6, self.c6, sign=+1)
        vr5T = self.xi_interface(self.ux4, self.ut4, self.c4, sign=+1)

        self.region_interfaces = [vr2H, vr2T, self.ux3, vr5T, vr5H]
        self.endInit()