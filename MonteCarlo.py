#!/home/arijitdutta/anaconda2/python2

import os
import sys
import random
import numpy as np

L = int(sys.argv[1])
sweep = int(sys.argv[2])
T = int(sys.argv[3])
J = 1.0
h = 0.0

out_config = './visualise/L' + str(L) + '/config_T' + str(T) + '_s' + str(sweep) + '.txt'
out_obs = './visualise/L' + str(L) + '/obs_T' + str(T) + '_s' + str(sweep) + '.txt'

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

def CalcLocalEnergyDiff(local_site, S_local_prev, S):
    deltaE = 0.0
    xcords = np.zeros(2)
    ycords = np.zeros(2)

    x = local_site / L;
    y = local_site % L;

    ## Periodic Boundary ##
    if (x > 0) and (x < L-1):
        xcords = [x-1, x+1]
    else:
        if (x == 0):
            xcords = [L-1, 1]
        if (x == L-1):
            xcords = [L-2, 0]

    if (y > 0) and (y < L-1):
        ycords = [y-1, y+1]
    else:
        if y == 0:
            ycords = [L-1, 1]
        if y == L-1:
            ycords = [L-2, 0]

    deltaE = (-J * (S[x, ycords[0]] + S[x, ycords[1]] + S[xcords[0], y] + S[xcords[1], y]) - h) * (S[x, y] - S_local_prev)
    return deltaE

def MonteCarlo(S):
    conf = open(out_config, 'w')
    obs = open(out_obs, 'w')

    for s in range(sweep):
        for i in range(L**2):
            ## Select a random site
            j = random.randint(0, (L**2)-1)
            x = j / L
            y = j % L
            ## Remember the current local configuration
            S_local_prev = S[x, y]
            ## Update the site
            S[x, y] = (-1)**random.randint(0, L**2)
            deltaE = CalcLocalEnergyDiff(j, S_local_prev, S)
            ## Metropolis decision tree
            if deltaE > 0.0:
                toss = random.random()
                beta = 1.0 / T
                if np.exp(-beta * deltaE) < toss:
                    S[x, y] = S_local_prev

        M = np.mean(S)
        E = CalculateEnergy(S)
        np.savetxt(conf, S.ravel(), fmt='%d', delimiter=',', newline=' ')
        np.savetxt(obs, np.c_[M, E], fmt='%f', newline=' ')
        conf.write('\n')
        obs.write('\n')

random.seed()
# Initialise a random spin configuration
S = np.array([[(-1)**random.randint(0, L**2) for i in range(L)] for j in range(L)])
MonteCarlo(S)
