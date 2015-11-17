import numpy as np

__all__ = ['Constraint', 'ConstPosition', 'Sinusoidal']

class Constraint(object):
    
    def __init__(self):
        self.time = 0
    
    def time(self, value):
        self.time = value        
    
    def check_constraint(self):
        pass                        #go to subclasses
    
class ConstPosition(Constraint):
    """Constrain body to a constant position"""
        
    def __init__(self, c):
        Constraint.__init__(self) #needed for time?
        self._constant = c
        
    def check_constraint(self):
        return (self._constant, 0)  

 
class Sinusoidal(Constraint):
        """Constrain body to a sinusoidal motion (e.g. pitching, heaving)"""
        
        def __init__(self, ampl, freq, phase):
            Constraint.__init__(self) #needed for time?
            self._ampl = ampl
            self._freq = freq
            self._phase = phase
    
        def check_constraint(self):
            x = self._ampl * np.sin(self._freq * self.time + self._phase)
            xdot = self._ampl * self._freq * np.cos(self._freq * self.time + self._phase)
            return (x, xdot)
        

    
    