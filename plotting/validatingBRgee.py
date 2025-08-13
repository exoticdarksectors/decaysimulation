import numpy as np
import math
from scipy.integrate import quad

# constants
mpi     = 0.1349768          # GeV
alpha   = 0.0072973526       # 1/137
PI      = math.pi
epsilon = 1.0
BrPi2gg = 1e9  # potential issue Pi2gg = 0.98823?

# differential branching ratio with BrPi2gg as input
def ddBrPi2gxx(s, theta, mchi, BrPi2gg):
    return (math.sin(theta)
            * epsilon*epsilon * alpha / (4.0*PI*s)
            * (1.0 - s/(mpi*mpi))**3
            * math.sqrt(1.0 - 4.0*mchi*mchi/s) * (2-(1 - 4.0*mchi*mchi/s)*math.sin(theta)**2) * BrPi2gg)

# Integrating function
def integrate_branching(mchi, BrPi2gg):
    s_min, s_max = 4*mchi*mchi, mpi*mpi

    # inner integral over θ for a fixed s
    def I_theta(s):
        return quad(lambda th: ddBrPi2gxx(s, th, mchi, BrPi2gg),
                    0.0, PI, epsabs=0, epsrel=1e-6)[0]

    # outer integral over s
    return quad(lambda s: I_theta(s),
                s_min, s_max, epsabs=0, epsrel=1e-6)[0]

# testing for electron branching ratio
if __name__ == "__main__":

    mchi = 0.000511 # GeV  (electron mass)

    BR_placeholder = integrate_branching(mchi, 1e9)
    BR_physical    = integrate_branching(mchi, 0.98823)

    print(f"mchi = {mchi} GeV")
    print(f"  BR  (placeholder) : {BR_placeholder:.4e}")
    print(f"  BR  (electron)    : {BR_physical:.4e}")
    print(f"  Ratio (placeholder/electron): {BR_placeholder/BR_physical:.3e}")