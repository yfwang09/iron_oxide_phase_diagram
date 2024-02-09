###########################################################################
# Calculate the thermodynamic equilibrium of H2O/H2 system in iron reduction
# process.
##############################################################################
import numpy as np

###########################################################################
# Constants
##############################################################################

# Gas constant
R = 8.314 # J/(mol K)

# Standard Gibbs free energy of Fe2O3
# (https://webbook.nist.gov/cgi/cbook.cgi?ID=C1317608&Mask=2)
def G_Fe2O3(T, Gibbs_only=True):
    H = -825.5032 # kJ/mol
    if 298 <= T and T < 950:
        A = 93.43834
        B = 108.3577
        C = -50.86447
        D = 25.58683
        E = -1.611330
        F = -863.2094
        G = 161.0719
    elif 950 <= T and T < 1050:
        A = 150.6240
        B = 0
        C = 0
        D = 0
        E = 0
        F = -875.6066
        G = 252.8814
    elif 1050 <= T and T < 2500:
        A = 110.9362
        B = 32.04714
        C = -9.192333
        D = 0.901506
        E = 0
        F = -843.1470
        G = 228.9738
    else:
        raise ValueError('Temperature out of range')
    t = T / 1000
    Hf = A*t + B*t**2/2 + C*t**3/3 + D*t**4/4 - E/t + F
    Sf = A*np.log(t) + B*t + C*t**2/2 + D*t**3/3 - E/(2*t**2) + G
    Gf = Hf - t*Sf

    if Gibbs_only:
        return Gf
    else:
        return Hf, Sf, Gf

# Standard Gibbs free energy of Fe3O4
# (https://webbook.nist.gov/cgi/formula?ID=C1309382&Mask=2)
def G_Fe3O4(T, Gibbs_only=True):
    H = -1120.89
    if 298 <= T and T < 900:
        A = 104.2096
        B = 178.5108
        C = 10.61510
        D = 1.132534
        E = -0.994202
        F = -1163.336
        G = 212.0585
    elif 900 <= T and T < 3000:
        A = 200.8320
        B = 1.586435e-7
        C = -6.661682e-8
        D = 9.452452e-9
        E = 3.186020e-8
        F = -1174.135
        G = 388.0790
    else:
        raise ValueError('Temperature out of range')
    t = T / 1000
    Hf = A*t + B*t**2/2 + C*t**3/3 + D*t**4/4 - E/t + F
    Sf = A*np.log(t) + B*t + C*t**2/2 + D*t**3/3 - E/(2*t**2) + G
    Gf = Hf - t*Sf

    if Gibbs_only:
        return Gf
    else:
        return Hf, Sf, Gf

# Standard Gibbs free energy of FeO
# (https://webbook.nist.gov/cgi/formula?ID=C1345251&Mask=2)
def G_FeO(T, phase='solid', Gibbs_only=True):
    if phase == 'solid' and 298 <= T and T < 1650:
        A = 45.75120
        B = 18.78553
        C = -5.952201
        D = 0.852779
        E = -0.081265
        F = -286.7429
        G = 110.3120
        H = -272.0441 # kJ/mol
    elif phase == 'liquid' and 1650 <= T and T < 5000:
        A = 68.19920
        B = -4.501232e-10
        C = 1.195227e-10
        D = -1.064302e-11
        E = -3.092680e-10
        F = -281.4326
        G = 137.8377
        H = -249.5321
    else:
        raise ValueError('Temperature out of range')
    t = T / 1000
    Hf = A*t + B*t**2/2 + C*t**3/3 + D*t**4/4 - E/t + F
    Sf = A*np.log(t) + B*t + C*t**2/2 + D*t**3/3 - E/(2*t**2) + G
    Gf = Hf - t*Sf

    if Gibbs_only:
        return Gf
    else:
        return Hf, Sf, Gf

# Standard Gibbs free energy of Fe
# (https://webbook.nist.gov/cgi/formula?ID=C7439896&Mask=2#Thermo-Condensed)
def G_Fe(T, phase='alpha', Gibbs_only=True):
    if phase == 'alpha' and 298 <= T and T < 1809:
        if 298 <= T and T < 700:
            A = 18.42868
            B = 24.64301
            C = -8.913720
            D = 9.664706
            E = -0.012643
            F = -6.573022
            G = 42.51488
            H = 0
        elif 700 <= T and T < 1042:
            A = -57767.65
            B = 137919.7
            C = -122773.2
            D = 38682.42
            E = 3993.080
            F = 24078.67
            G = -87364.01
            H = 0
        elif 1042 <= T and T < 1100:
            A = -325.8859
            B = 28.92876
            C = 0
            D = 0
            E = 411.9629
            F = 745.8231
            G = 241.8766
            H = 0
        elif 1100 <= T and T < 1809:
            A = -776.7387
            B = 919.4005
            C = -383.7184
            D = 57.08148
            E = 242.1369
            F = 697.6234
            G = -558.3674
            H = 0
    elif phase == 'delta' and 298 <= T and T < 1809:
        A = 23.97449
        B = 8.367750
        C = 0.000277
        D = -0.000086
        E = -0.000005
        F = 0.268027
        G = 62.06336
        H = 7.788015
    elif phase == 'liquid' and 1809 <= T and T < 3133.345:
        A = 46.02400
        B = -1.884667e-8
        C = 6.094750e-9
        D = -6.640301e-10
        E = -8.246121e-9
        F = -10.80543
        G = 72.54094
        H = 12.39502
    else:
        raise ValueError('Temperature out of range')
    t = T / 1000
    Hf = A*t + B*t**2/2 + C*t**3/3 + D*t**4/4 - E/t + F
    Sf = A*np.log(t) + B*t + C*t**2/2 + D*t**3/3 - E/(2*t**2) + G
    Gf = Hf - t*Sf

    if Gibbs_only:
        return Gf
    else:
        return Hf, Sf, Gf

# Standard Gibbs free energy of H2O
# (https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Mask=2)
def G_H2O(T, phase='gas', Gibbs_only=True):
    if phase == 'gas' and 300 <= T and T < 6000:
        H = -241.8264
        if 300 <= T and T < 1700:
            A = 30.09200
            B = 6.832514
            C = 6.793435
            D = -2.534480
            E = 0.082139
            F = -250.8810
            G = 223.3967
        elif 1700 <= T and T < 6000:
            A = 41.96426
            B = 8.622053
            C = -1.499780
            D = 0.098119
            E = -11.15764
            F = -272.1797
            G = 219.7809
    elif phase == 'liquid' and 298 <= T and T < 500:
        A = -203.6060
        B = 1523.290
        C = -3196.413
        D = 2474.455
        E = 3.855326
        F = -256.5478
        G = -488.7163
        H = -285.8307
    else:
        raise ValueError('Temperature out of range')
    t = T / 1000
    Hf = A*t + B*t**2/2 + C*t**3/3 + D*t**4/4 - E/t + F
    Sf = A*np.log(t) + B*t + C*t**2/2 + D*t**3/3 - E/(2*t**2) + G
    Gf = Hf - t*Sf

    if Gibbs_only:
        return Gf
    else:
        return Hf, Sf, Gf

# Standard Gibbs free energy of H2
# (https://webbook.nist.gov/cgi/cbook.cgi?ID=C1333740&Mask=2)
def G_H2(T, phase='gas', Gibbs_only=True):
    if phase == 'gas' and 298 <= T and T < 1000:
        A = 33.066178
        B = -11.363417
        C = 11.432816
        D = -2.772874
        E = -0.158558
        F = -9.980797
        G = 172.707974
    elif phase == 'gas' and 1000 <= T and T < 2500:
        A = 18.563083
        B = 12.257357
        C = -2.859786
        D = 0.268238
        E = 1.977990
        F = -1.147438
        G = 156.288133
    elif phase == 'gas' and 2500 <= T and T < 6000:
        A = 43.413560
        B = -4.293079
        C = 1.272428
        D = -0.096876
        E = -20.533862
        F = -38.515158
        G = 162.081354
    else:
        raise ValueError('Temperature out of range')
    t = T / 1000
    Hf = A*t + B*t**2/2 + C*t**3/3 + D*t**4/4 - E/t + F
    Sf = A*np.log(t) + B*t + C*t**2/2 + D*t**3/3 - E/(2*t**2) + G
    Gf = Hf - t*Sf

    if Gibbs_only:
        return Gf
    else:
        return Hf, Sf, Gf
