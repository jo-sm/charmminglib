"""
DOCME
"""
# fcp
# 10/27/2010


from numpy import log
from pychm.const.units import BOLTZMANN, AVOGADRO, JOULE2CAL, CAL2JOULE
from pychm.tools import Property


class DelGError(Exception):
    """
    Exception to raise when errors occur involving the DelG class.
    """
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


class DelG(object):
    """
    docstring for DelG
    """
    def __init__(self, timeSeries, temp):
        self.timeSeries = timeSeries
        self.stateDef = {}
        self.stateCount = {}
        self.temp = temp

    @Property
    def temp():
        doc = "The temp property."
        def fget(self):
            return self._temp
        def fset(self, value):
            value = float(value)
            if value <= 0.:
                raise DelGError('temp: temp must be specified in K,\
                                and greater than 0')
            else:
                self._temp = value
        return locals()

    def addState(self, name, less, greater):
        if less > greater:
            less, greater = greater, less
        self.stateDef[name] = (less, greater)

    def count(self):
        # Count state populations
        for key in self.stateDef.iterkeys():
            self.stateCount[key] = 0
        for frame in self.timeSeries:
            for key, value in self.stateDef.iteritems():
                if value[0] <= frame <= value[1]:
                    self.stateCount[key] += 1
                    break

    def get_DelG(self, state0, state1):
        if not self.temp:
            raise DelGError('get_DelG: temp not specified')
        BETA = 1. / (BOLTZMANN * self.temp)
        count0 = float(self.stateCount[state0])
        count1 = float(self.stateCount[state1])
        energy =  -1 * log(count0 / count1) / BETA  # Joule
        energy = energy / 1000.                 # Kilo Joule
        energy = energy * AVOGADRO              # Joule * Mole**-1
        energy = energy * JOULE2CAL             # KCal / mol
        return energy
