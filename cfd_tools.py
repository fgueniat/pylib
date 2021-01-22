'''
library to calculate nut, k, epsilon and omega.
Florimond Gueniat 2020
'''
import numpy as np

def estimate_dt(uinf,dx):
    return dx/uinf

def turbulence_intensity(uinf,Lpipe,nu):
    '''
    Turbulence intensity I
    https://www.cfd-online.com/Wiki/Turbulence_intensity
    '''
    Repipe =uinf*Lpipe/nu
    I = 0.16 * (Repipe**(-1/8.))
    return I

def turbulence_scale(Lpipe):
    '''
    Turbulencec scale
    https://www.cfd-online.com/Wiki/Turbulent_length_scale
    '''
    l=0.038*Lpipe
    return l

def nut(uinf,I=False,l=False,Lpipe=False,nu=False):
    '''
    binding to nu tilda
    Needs the turbulence intensity I.
    Altenatively, I is calculated using the wind tunnel diameter Lpipe and the kinematic viscosity nu
    Needs the turbulent scale l.
    Alternatively, l is calculated using the wind tunnel diameter Lpipe

    https://www.cfd-online.com/Wiki/Turbulence_free-stream_boundary_conditions
    '''
    return nu_tilda(uinf,I,l,Lpipe,nu)

def nu_tilda(uinf,I=False,l=False,Lpipe=False,nu=False):
    '''
    nu tilda
    Needs the turbulence intensity I.
    Altenatively, I is calculated using the wind tunnel diameter Lpipe and the kinematic viscosity nu
    Needs the turbulent scale l.
    Alternatively, l is calculated using the wind tunnel diameter Lpipe

    https://www.cfd-online.com/Wiki/Turbulence_free-stream_boundary_conditions
    '''
    if I is False:
        if Lpipe is False:
            print('please provide the diameter of the wind tunnel')
            return -1
        if nu is False:
            print('please provide the kinematic viscosity')
            return -1
        I = turbulence_intensity(uinf,Lpipe,nu)
    if l is False:
        if Lpipe is False:
            print('please provide the diameter of the wind tunnel')
        l = turbulence_scale(Lpipe)

    nut = np.sqrt(3/2) * uinf * I * l
    return nut


def kte(uinf,I=False,Lpipe=False,nu=False):
    '''
    binding to turbulent energy k
    Needs the turbulence intensity I.
    Alternatively, I is calculated using the wind tunnel diameter Lpipe and the kinematic viscosity nu

    https://www.cfd-online.com/Wiki/Turbulence_free-stream_boundary_conditions
    '''
    return turbulent_energy(uinf,I,Lpipe,nu)

def turbulent_energy(uinf,I=False,Lpipe=False,nu=False):
    '''
    turbulent energy k
    Needs the turbulence intensity I.
    Alternatively, I is calculated using the wind tunnel diameter Lpipe and the kinematic viscosity nu

    https://www.cfd-online.com/Wiki/Turbulence_free-stream_boundary_conditions
    '''
    if I is False:
        if Lpipe is False:
            print('please provide the diameter of the wind tunnel')
            return -1
        if nu is False:
            print('please provide the kinematic viscosity')
            return -1
        I = turbulence_intensity(uinf,Lpipe,nu)

    k = (3./2) * ((uinf*I)**2)
    return k

def epsilon(uinf,k=False,l=False,I=False,Lpipe=False,nu=False):
    '''
    binding to dissipation rate epsilon
    Needs the turbulence energy k.
    Alternatively, k is calculated using the turbulent intensity I.
    As a second alternatively, I can calculated using the wind tunnel diameter Lpipe and the kinematic viscosity nu
    Needs the turbulent scale l.
    Alternatively, l is calculated using the wind tunnel diameter Lpipe
    
    https://www.cfd-online.com/Wiki/Turbulence_free-stream_boundary_conditions
    '''
    return dissipation_rate(uinf,k,I,l,Lpipe,nu)

def dissipation_rate(uinf,k=False,I=False,l=False,Lpipe=False,nu=False):
    '''
    dissipation rate epsilon
    Needs the turbulence energy k.
    Alternatively, k is calculated using the turbulent intensity I.
    As a second alternatively, I can calculated using the wind tunnel diameter Lpipe and the kinematic viscosity nu
    Needs the turbulent scale l.
    Alternatively, l is calculated using the wind tunnel diameter Lpipe

    https://www.cfd-online.com/Wiki/Turbulence_free-stream_boundary_conditions
    '''
    if k is False:
        if I is False:
            if Lpipe is False:
                print('please provide the diameter of the wind tunnel')
                return -1
            if nu is False:
                print('please provide the kinematic viscosity')
                return -1
            I = turbulence_intensity(uinf,Lpipe,nu)
        k = turbulent_energy(uinf,I)
    if l is False:
        if Lpipe is False:
            print('please provide the diameter of the wind tunnel')
        l = turbulence_scale(Lpipe)

    epsilon = cmu * (k**(3./2)) / l
    return epsilon


def omega(uinf,k=False,l=False,I=False,Lpipe=False,nu=False):
    '''
    binding to specific dissipation rate omega
    Needs the turbulence energy k.
    Alternatively, k is calculated using the turbulent intensity I.
    As a second alternatively, I can calculated using the wind tunnel diameter Lpipe and the kinematic viscosity nu
    Needs the turbulent scale l.
    Alternatively, l is calculated using the wind tunnel diameter Lpipe
    
    https://www.cfd-online.com/Wiki/Turbulence_free-stream_boundary_conditions
    '''
    return specific_dissipation_rate(uinf,k,I,l,Lpipe,nu)

def specific_dissipation_rate(uinf,k=False,I=False,l=False,Lpipe=False,nu=False):
    '''
    dissipation rate omega
    Needs the turbulence energy k.
    Alternatively, k is calculated using the turbulent intensity I.
    As a second alternatively, I can calculated using the wind tunnel diameter Lpipe and the kinematic viscosity nu
    Needs the turbulent scale l.
    Alternatively, l is calculated using the wind tunnel diameter Lpipe

    https://www.cfd-online.com/Wiki/Turbulence_free-stream_boundary_conditions
    '''
    if k is False:
        if I is False:
            if Lpipe is False:
                print('please provide the diameter of the wind tunnel')
                return -1
            if nu is False:
                print('please provide the kinematic viscosity')
                return -1
            I = turbulence_intensity(uinf,Lpipe,nu)
        k = turbulent_energy(uinf,I)
    if l is False:
        if Lpipe is False:
            print('please provide the diameter of the wind tunnel')
        l = turbulence_scale(Lpipe)

    omega = np.sqrt(k) / l
    return omega




