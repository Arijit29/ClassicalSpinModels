#!/usr/bin/python2.6

import sys
import numpy as np
import matplotlib.pyplot as plt

L = int(sys.argv[1])
sweep = int(sys.argv[2])
T = int(sys.argv[3])

plt.ion()
fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(16, 6))

in_file = './visualise/L' + str(L) + '/obs_T' + str(T) + '_s' + str(sweep) + '.txt'

# Assuming sweep >= 100000
Bins = np.array([1, 10, 20, 40, 50, 80, 100, 125, 200, 250, 400, 500, 800, 1000])
Nbins = np.array([sweep / b for b in Bins])

M_av = 0
E_av = 0
M_std = np.zeros(len(Bins))
E_std = np.zeros(len(Bins))

print 'Reading configurations from file -> ', in_file
obs = np.zeros(sweep)
M_config = np.zeros(sweep)
E_config = np.zeros(sweep)
M_bins = [np.zeros(Nbins[i]) for i in range(len(Nbins))]
E_bins = [np.zeros(Nbins[i]) for i in range(len(Nbins))]

s = 0
cbin = np.zeros(len(Nbins), dtype=int)
ax0.set_ylim(-1, 1)
ax0.set_yticks(np.arange(-1, 1, 0.2))
for line in open(in_file).readlines():
    temp = np.fromstring(line, dtype=float, sep=' ')
    M_config[s] = temp[0]
    E_config[s] = temp[1]
    M_av += M_config[s]
    E_av += E_config[s]
    if s <= 500:
        ax0.scatter(s, M_av/(s+1), marker='o', label='<M>', color='r')
        ax0.scatter(s, E_av/(s+1), marker='s', label='<E>', color='b')
        plt.pause(0.0005)
    if s == 0:
        ax0.grid()
        ax0.legend()
        ax0.set_xlabel('MC time')
    for i in range(len(Nbins)):
    	M_bins[i][cbin[i]] += M_config[s]
    	E_bins[i][cbin[i]] += E_config[s]
    s = s + 1
    for i in range(len(Bins)):
    	if (s % Bins[i]) == 0:
            cbin[i] = cbin[i] + 1

for i in range(len(Nbins)):
    M_bins[i] /= Bins[i]
    E_bins[i] /= Bins[i]
    for j in range(Nbins[i]):
        M_std[i] += (M_bins[i][j])**2
        E_std[i] += (E_bins[i][j])**2
    M_std[i] /= Nbins[i]
    E_std[i] /= Nbins[i]
    M_std[i] -= (np.mean(M_bins[i]))**2
    E_std[i] -= (np.mean(E_bins[i]))**2
    M_std[i] = (M_std[i]/ (Nbins[i]-1))**(0.5)
    E_std[i] = (E_std[i]/ (Nbins[i]-1))**(0.5)

out_file = './visualise/L' + str(L) + '/std_T' + str(T) + '_s' + str(sweep) + '.txt'
np.savetxt(out_file, np.c_[Bins, M_std, E_std], delimiter='\t')
ax0.plot(E_av, color='b')
ax1.plot(Bins, M_std, marker='o', label='$\sigma^{M}$', color='r')
ax1.plot(Bins, E_std, marker='s', label='$\sigma^{E}$', color='b')
ax1.set_xlabel('Bin size')
# ax.set_ylabel('$\sigma$')
ax1.set_title('T=' + str(T) + ', sweep=' + str(sweep))
ax1.grid()
ax1.legend()
ax1.set_ylim(0, 0.001)
plt.show()

while True:
    plt.pause(0.05)
