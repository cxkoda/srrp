from .ContactDiscontinuity import ContactDiscontinuity
import numpy as np
from .Rarefaction import Rarefaction
from .Shock import Shock
from .State import State


class Wavefan:
    def __init__(self, states, waves, reversed=False):
        assert (len(states) == len(waves) + 1)
        self.states = states
        self.waves = waves
        self.regionBoundaries, self.regionStates = Wavefan.__build(states, waves)

    @staticmethod
    def __build(states, waves):
        regionBoundaries = [-np.inf]
        regionStates = []

        regionStates.append(lambda xi: states[0])
        for idx, _ in enumerate(waves):
            Wavefan.__appendWave(regionBoundaries, regionStates, waves[idx])
            regionStates.append(lambda xi, idx=idx: states[idx + 1])

        regionBoundaries.append(np.inf)
        return regionBoundaries, regionStates

    @staticmethod
    def __appendWave(regionBoundaries, regionStates, wave):
        if isinstance(wave, Shock):
            regionBoundaries.append(wave.speed)
        elif isinstance(wave, ContactDiscontinuity):
            regionBoundaries.append(wave.speed)
        elif isinstance(wave, Rarefaction):
            regionBoundaries.append(min(wave.speedTail, wave.speedHead))
            regionStates.append(wave.computeRarefactionState)
            regionBoundaries.append(max(wave.speedTail, wave.speedHead))
        else:
            raise RuntimeError("Unknown wave")

    @staticmethod
    def __getRegionCondlist(xi, regionBoundaries):
        condlist = []
        condlist.append(xi <= regionBoundaries[0])
        for left, right in zip(regionBoundaries[:-1], regionBoundaries[1:]):
            condlist.append((left < xi) * (xi <= right))
        condlist.append(regionBoundaries[-1] < xi)
        return condlist

    # @staticmethod
    # def __getRegionIndex(xi, regionBoundaries):
    #     condlist = []
    #     condlist.append(xi <= regionBoundaries[0])
    #     for left, right in zip(regionBoundaries[:-1], regionBoundaries[1:]):
    #         condlist.append((left < xi) * (xi <= right))
    #     condlist.append(regionBoundaries[-1] < xi)
    #     return condlist

    def getRegIdx(self, x, t):
        xi = self.get_xi(x, t)

        for idx, interface in enumerate(self.region_interfaces):
            if interface is not None and xi <= interface:
                return idx + 1

        return 6

    def getRegionIndex(self, xi):
        return np.digitize(xi, self.regionBoundaries) - 1


    def __getState(self, xi):
        idx = self.getRegionIndex(xi)
        return self.regionStates[idx](xi)

    def getState(self, xis):
        doGet = np.vectorize(self.__getState)
        states = doGet(xis)

        try:
            statesAggregated = type(states[0])()
            for var in vars(statesAggregated):
                setattr(statesAggregated, var, np.array([getattr(state, var) for state in states]))
            return statesAggregated
        except:
            return states


    def __str__(self):
        tmp = [str(self.states[0])]
        for wave, state in zip(self.waves, self.states):
            tmp.append(str(wave))
            tmp.append(str(state))
        return '\n'.join(tmp)

    def __repr__(self):
        return str(self)