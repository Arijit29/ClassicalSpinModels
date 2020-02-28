#!/usr/bin/python2.6

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

L = int(sys.argv[1])
sweep = int(sys.argv[2])
Nconf = int(sys.argv[3])
J = 1.0
h = 0.0
Temperature = [5, 4.5, 4, 3.75, 3.5, 3.25, 3, 2.8, 2.6, 2.5, 2.4, 2.3, 2.2, 2.1, 2, 1.8, 1.6, 1.4, 1 , 0.01]

out_dir = './L' + str(L) + '/'

M_av = np.zeros_like(Temperature)
E_av = np.zeros_like(Temperature)
MagSusp = np.zeros_like(Temperature)
Cv = np.zeros_like(Temperature)

def CalculateEnergy(S):
    #H = np.zeros(shape=(L, L))
    Energy = 0.0
    ## Bulk Hamiltonian ##
    for i in range(1, L-1):
        for j in range(1, L-1):
            Energy += (-J * S[i, j] * (S[i, j-1] + S[i, j+1] + S[i-1, j] + S[i+1, j])) + (-h * S[i, j]);

    ## Periodic Boundary Conditions ##
    for j in range(1, L-1):
        Energy += (-J * S[0, j] * (S[0, j-1] + S[0, j+1] + S[L-1, j] + S[1, j])) + (-h * S[0, j]);
        Energy += (-J * S[j, 0] * (S[j-1, 0] + S[j+1, 0] + S[j, L-1] + S[j, 1])) + (-h * S[j, 0]);
        Energy += (-J * S[j, L-1] * (S[j-1, L-1] + S[j+1, L-1] + S[j, L-2] + S[j, 0])) + (-h * S[j, L-1]);
        Energy += (-J * S[L-1, j] * (S[L-1, j-1] + S[L-1, j+1] + S[L-2, j] + S[0, j])) + (-h * S[L-1, j]);

    Energy += (-J * S[0, 0] * (S[0, L-1] + S[0, 1] + S[L-1, 0] + S[1, 0])) + (-h * S[0, 0]);
    Energy += (-J * S[L-1, 0] * (S[L-1, L-1] + S[L-1, 1] + S[L-2, 0] + S[0, 0])) + (-h * S[L-1, 0]);
    Energy += (-J * S[L-1, L-1] * (S[L-1, L-2] + S[L-1, 0] + S[L-2, L-1] + S[0, L-1])) + (-h * S[L-1, L-1]);
    Energy += (-J * S[0, L-1] * (S[0, L-2] + S[0, 0] + S[L-1, L-1] + S[1, L-1])) + (-h * S[0, L-1]);
    #*********************************#
    return (Energy / float(L*L))

b = 0
for T in Temperature:
    filename = out_dir + 'config_T' + str(T) + '.txt'
    print 'Reading configurations from file -> ', filename
    config = np.zeros(shape=(Nconf, L*L))
    S = np.zeros(shape=(L, L))
    M_config = np.zeros(Nconf)
    E_config = np.zeros(Nconf)
    s = 0
    cnt = 0
    save = (sweep / (2 * Nconf))
    for line in reversed(open(filename).readlines()):
        if (s < sweep/2):
            if (s % save == 0) and (cnt < Nconf):
                config[cnt, :] = np.fromstring(line, dtype=float, sep='\t')
                for i in range(0, L**2, 1):
                    S[i/L, i%L] = config[cnt, i]

                M_config[cnt] = np.mean(config[cnt,:])
                E_config[cnt] = CalculateEnergy(S)
                cnt = cnt + 1
        else:
            break
        s = s + 1

    M_av[b] = abs(np.mean(M_config))
    E_av[b] = np.mean(E_config)
    MagSusp[b] = (np.var(M_config) / T)
    Cv[b] = (np.var(E_config)/ T**2)
    b = b + 1

np.savetxt(out_dir+'obs.txt', np.c_[Temperature, M_av, E_av, MagSusp, Cv], fmt='%.3E', delimiter='\t', newline=os.linesep)

fig, ((ax1, ax2), (bx1, bx2)) = plt.subplots(2, 2, sharex=True)
ax1.set_xlim(0, 5)
ax2.set_xlim(0, 5)
bx1.set_xlim(0, 5)
bx2.set_xlim(0, 5)

ax1.plot(Temperature, M_av, marker='o')
ax2.plot(Temperature, MagSusp, marker='^')
bx1.plot(Temperature, E_av, marker='o')
bx2.plot(Temperature, Cv, marker='^')

plt.show()

