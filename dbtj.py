
import numpy as _np
import numba as _nb
from scipy.constants import electron_volt as _e, Boltzmann as _kB

@_nb.njit
def gamma(dE, R, T=4.2, val=1e-1):
    # type: (numpy.ndarray, float, float, float) -> numpy.ndarray
    '''

    Fermi's golden rule tunneling rate.

    Parameters
    ----------
    dE: numpy.ndarray
        Energy change in Joules.
    R: float
        Gap resistance in Ohms.
    T: float, default is 4.2
        System temperature in Kelvin.
    val: float, default is 1e-1
        Value to replace unrealistic tunneling rates with to avoid ZeroDivisionError.
        Keep small.
    
    Returns
    -------
    gamma: numpy.ndarray
        Tunneling rates.
    
    '''
    # format input
    dE = _np.asarray(dE, dtype=_np.float_)
    # calculate model
    out = (1./(R*_e**2.)) * -dE/(1.-_np.exp(dE/(_kB*T)))
    # substitute in val is dE>0 -> gamma would be small or zero
    # avoids ZeroDivisionError
    out[dE>0.] = val
    return out

@_nb.njit
def dbtj(bias, R1, R2, C1, C2, Q0, T=4.2, nmax=20, val=1e-1):
    # type: (numpy.ndarray, float, float, float, float, float, float, int, float) -> numpy.ndarray
    '''
    
    Double Barrier Tunneling Junction Model.

    [1] Hanna, Tinkham (1991), Phys. Rev. B 44 (11), p. 5919-5922, DOI:10.1103/PhysRevB.44.5919
    
    Parameters
    ----------
    bias: array-like
        Bias values to simulate.
    R1, R2: float
        Resistance of junctions 1 and 2 in Ohms.
    C1, C2: float
        Capacitance of junctions 1 and 2 in Farads.
    Q0: float, value lies within (-0.5<=Q0<0.5)
        Fractional charge on center electrode in units of elementary charge.
    T: float, default is 4.2
        System temperature in Kelvin
    nmax: int, default is 20
        Number of electrons on center electrode to consider.
    val: float, default is 1e-1
        Value to replace unrealistic tunneling rates with to avoid ZeroDivisionError.
        Keep small.
    
    Returns
    -------
    I: numpy.ndarray
        Calculated tunneling current.
    
    '''
    # format input
    bias = _np.asarray(bias, dtype=_np.float_)
    # check valid Q0
    assert (-1./2)<=Q0<(1./2), 'Q0 must be in range (-0.5<=Q0<0.5).'
    
    # total capacitance
    Csum = C1 + C2
    
    # number of electrons on center electrode to consider
    n = _np.arange(-nmax, nmax+1)
    
    # current array holders
    I1 = _np.zeros_like(bias, dtype=_np.float_)
    # uncomment below to return both I1 and I2
    # I2 = _np.zeros_like(I1)

    # for each bias
    for i, V in enumerate(bias):
        # calculate ±(dE1, dE2) according to equations
        dE1_plus = _e/Csum * (_e/2. + _e*(n-Q0) + C2*V)
        dE1_minus = _e/Csum * (_e/2. - _e*(n-Q0) - C2*V)

        dE2_plus = _e/Csum * (_e/2. + _e*(n-Q0) - C1*V)
        dE2_minus = _e/Csum * (_e/2. - _e*(n-Q0) + C1*V)

        # calculate tunnelling rate accoring to Fermi's golden rule approx.
        gamma1_plus = gamma(dE1_plus, R1, T, val)
        gamma1_minus = gamma(dE1_minus, R1, T, val)

        gamma2_plus = gamma(dE2_plus, R2, T, val)
        gamma2_minus = gamma(dE2_minus, R2, T, val)

        # make sure that gamma is nonzero
        gamma1_plus[dE1_plus>0.] = val
        gamma1_minus[dE1_minus>0.] = val
        gamma2_plus[dE2_plus>0.] = val
        gamma2_minus[dE2_minus>0.] = val

        # calculate sigma and normalise -> sigma(n)/sigma(n+1)
        sigma = _np.ones_like(n, dtype=_np.float_)
        # sigma(n+1) is dependent on sigma(n)
        for j in range(sigma[1:].size):
            sigma[j+1] = sigma[j] * (gamma1_plus[j] + gamma2_plus[j])\
                                  / (gamma1_minus[j+1] + gamma2_minus[j+1])
        # normalise sigma
        sigma /= sigma.sum()
        
        # calculate DBTJ current
        I1[i] = _e * _np.sum(sigma * (gamma1_minus-gamma1_plus))
        # I2[i] = _e * _np.sum(sigma * (gamma2_plus-gamma2_minus))
    
    out = I1
    # uncomment below to return both I1 and I2
    # out = (I1, I2)
    return out
