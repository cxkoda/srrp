from .SRHDRiemannSolver import SRHDRiemannSolver

class SRHDRiemannSolverCronosPy:
    '''
    Interface class for CronosPy
    '''

    def set_pars(self, _rhoL, _pL, _vxL, _vtL, _rhoR, _pR, _vxR, _vtR, _gamma, _time):
        self.time = _time
        # Attention! Argument order is different!
        self.solution = SRHD_RiemannSolver(_rhoL, _vxL, _vtL, _pL, _rhoR, _vxR, _vtR, _pR, _gamma, verbose=True)

    def reset_time(self, _time):
        self.time = _time

    def set_field(self, _fieldName):
        self.fieldName = _fieldName

    def set_mod(self, _modName):
        self.modName = _modName

    def get_yvals(self, xData):
        print(" Name: " + self.fieldName)

        vals = 0

        if self.fieldName == "Density":
                vals = [self.get_rho(x) for x in xData]
        elif self.fieldName == "Velocity":
                vals = [self.get_vel(x) for x in xData]
        elif self.fieldName == "Etherm":
                vals = [self.get_Etherm(x) for x in xData]
        else:
                vals = [self.get_rho(x) for x in xData]

        vals = np.array(vals)

        if self.modName == ' Abs ':
            if len(vals.shape) > 1:
                vals = np.linalg.norm(vals, axis=1)
            else:
                vals = np.abs(vals)
        elif self.modName == ' v_x ':
            vals = vals[:, 0]
        elif self.modName == ' v_y ':
            vals = vals[:, 1]

        return vals

    def get_rho(self, xPos):
        rho, ux, ut, pres = self.get_allValues(xPos)

        return rho

    def get_vel(self, xPos):
        rho, ux, ut, pres = self.get_allValues(xPos)

        return [ux, ut]

    def get_pres(self, xPos):
        rho, ux, ut, pres = self.get_allValues(xPos)

        return pres

    def get_Etherm(self, xPos):
        rho, ux, ut, pres = self.get_allValues(xPos)

        return pres / (self.gamma - 1)

    def get_Temp(self, xPos):
        rho, ux, ut, pres = self.get_allValues(xPos)

        return pres / rho

    def get_allValues(self, xPos):
        return self.solution(xPos, self.time)
