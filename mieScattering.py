import numpy as np
def mie_scattering(r, f, epsp, epsb):
    """
    Calculates the Mie scattering efficiency factors resulting from spherical
    inclusions embedded in a medium.
    
    Parameters
    ----------
    r : float
        Particle radius (m)
    f : float
        Frequency (Hz)
    epsp : complex
        Particle relative permittivity
    epsb : complex
        Background relative permittivity
    
    Returns
    -------
    Es : float
        Scattering efficiency factor
    Ea : float
        Absorption efficiency factor
    Ee : float
        Extinction efficiency
    Eb : float
        Backscattering efficiency
    
    Source
    ------
    Ulaby and Long (2014)
    
    Author
    ------
    Ulaby and Long (2014) modified by Natalie Wolfenbarger
    nswolfen@gmail.com
    Converted to Python
    """
    
    c = 3e8  # m/s
    epsb_real = np.real(epsb)
    
    np_refract = np.sqrt(epsp)  # index of refraction of spherical particle
    nb = np.sqrt(epsb)  # index of refraction of background medium
    
    n = np_refract / nb  # relative index of refraction (8.31a)
    
    lambda_wave = c / f
    chi = (2 * np.pi * r / lambda_wave) * np.sqrt(epsb_real)  # normalized circumference in background (8.31b)
    
    # ========== Calculate Es (Scattering efficiency) ==========
    # Values of W0 and W-1
    W_1 = np.sin(chi) + 1j * np.cos(chi)  # (8.35a)
    W_2 = np.cos(chi) - 1j * np.sin(chi)  # (8.35b)
    
    # Value of A0
    A = 1.0 / np.tan(n * chi)  # (8.37)
    
    oldSum = 0
    pdiff = 1
    l = 1
    
    while pdiff >= 0.001:
        W = (2 * l - 1) / chi * W_1 - W_2  # (8.34)
        
        A = -l / (n * chi) + (l / (n * chi) - A) ** (-1)  # (8.36)
        
        a = ((A / n + l / chi) * np.real(W) - np.real(W_1)) / ((A / n + l / chi) * W - W_1)  # (8.33a)
        b = ((n * A + l / chi) * np.real(W) - np.real(W_1)) / ((n * A + l / chi) * W - W_1)  # (8.33b)
        
        sigma = (2 * l + 1) * (np.abs(a) ** 2 + np.abs(b) ** 2)  # (8.32a)
        newSum = oldSum + sigma
        
        l = l + 1
        W_2 = W_1
        W_1 = W
        
        pdiff = np.abs((newSum - oldSum) / newSum) * 100
        oldSum = newSum
    
    Es = 2 / (chi) ** 2 * newSum  # (8.32a)
    
    # ========== Calculate Ee (Extinction efficiency) ==========
    # Values of W0 and W-1
    W_1 = np.sin(chi) + 1j * np.cos(chi)  # (8.35a)
    W_2 = np.cos(chi) - 1j * np.sin(chi)  # (8.35b)
    
    # Value of A0
    A = 1.0 / np.tan(n * chi)  # (8.37)
    
    oldSum = 0
    pdiff = 1
    l = 1
    
    while pdiff >= 0.001:
        W = (2 * l - 1) / chi * W_1 - W_2  # (8.34)
        
        A = -l / (n * chi) + (l / (n * chi) - A) ** (-1)  # (8.36)
        
        a = ((A / n + l / chi) * np.real(W) - np.real(W_1)) / ((A / n + l / chi) * W - W_1)  # (8.33a)
        b = ((n * A + l / chi) * np.real(W) - np.real(W_1)) / ((n * A + l / chi) * W - W_1)  # (8.33b)
        
        sigma = (2 * l + 1) * np.real(a + b)  # (8.32b)
        newSum = oldSum + sigma
        
        l = l + 1
        W_2 = W_1
        W_1 = W
        
        pdiff = np.abs((newSum - oldSum) / newSum) * 100
        oldSum = newSum
    
    Ee = 2 / (chi) ** 2 * newSum  # (8.32b)
    
    # ========== Calculate Eb (Backscattering efficiency) ==========
    # Values of W0 and W-1
    W_1 = np.sin(chi) + 1j * np.cos(chi)  # (8.35a)
    W_2 = np.cos(chi) - 1j * np.sin(chi)  # (8.35b)
    
    # Value of A0
    A = 1.0 / np.tan(n * chi)  # (8.37)
    
    oldSum = 0
    pdiff = 1
    l = 1
    
    while pdiff >= 0.001:
        W = (2 * l - 1) / chi * W_1 - W_2  # (8.34)
        
        A = -l / (n * chi) + (l / (n * chi) - A) ** (-1)  # (8.36)
        
        a = ((A / n + l / chi) * np.real(W) - np.real(W_1)) / ((A / n + l / chi) * W - W_1)  # (8.33a)
        b = ((n * A + l / chi) * np.real(W) - np.real(W_1)) / ((n * A + l / chi) * W - W_1)  # (8.33b)
        
        sigma = (-1) ** l * (2 * l + 1) * (a - b)  # (8.32b)
        newSum = oldSum + sigma
        
        l = l + 1
        W_2 = W_1
        W_1 = W
        
        pdiff = np.abs((newSum - oldSum) / newSum) * 100
        oldSum = newSum
    
    Eb = 1 / (chi) ** 2 * np.abs(newSum) ** 2  # (8.40)
    
    Ea = Ee - Es  # (8.28b)
    
    return Es, Ea, Ee, Eb