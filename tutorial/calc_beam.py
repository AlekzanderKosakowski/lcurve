import numpy as np

def Beaming_factor(T,l): # Temp, wavelength(A)
    """
    Calculate the relativisitic beaming factor using approximation by Loeb and Gaudi

    input:
    T: float
    temperature in K
    l: float
    wavelength in AA

    """

    # as a test; T=5700, l=6000AA should give B=4.3
    l = l*10**-10 # convert to meters
    e = 0.1*10**-9
    B = 5 + (np.log(Fl_planck(l+e,T))-np.log(Fl_planck(l,T)) ) / (np.log(l+e)-np.log(l))

    return B

def Fl_planck(l,T):

    '''
    input :
    l : float or array-like
    wavelength in meters
    T : float or array-like
    temperature in Kelvin

    output : float or array-like
    the effective radius divided by the semimajor axis
    '''

    h = 6.626070040*10**-34
    kb = 1.38064852*10**-23
    c = 299792458.

    return 2.*h*l**-5.*c**-2. * (np.exp((h*c)/(l*kb*T))-1)**-1


if __name__ == "__main__":

  import sys

  t = int(sys.argv[1])
  filter = sys.argv[2]

  filters = {'g': 4470,
             'r': 6231,
             'i': 7625,}


  l = filters[filter]
#  l = 4470
#  t = 17500
#  l = 4470

  print(Beaming_factor(t,l))
