import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['text.usetex'] = True

ranL = [24, 36, 48]
rangeT = [5, 4.5, 4.0, 3.75, 3.5, 3.25, 3.0, 2.8, 2.6, 2.5, 2.4, 2.3, 2.2, 2.1, 2.0, 1.8, 1.6, 1.4, 1.0 , 0.01]
k0 = int(0)
markers = ['o', 's', '^']
colors = ['tab:green', 'tab:blue', 'tab:red']
fig, ax = plt.subplots(figsize=(4.5,4.5))

for l in range(len(ranL)):
	L = ranL[l]
	path = 'L'+str(L)+'/'

	SkPeak = np.array([]).reshape(0,2)
	for t in range(len(rangeT)):
		T = rangeT[t]
		if T == int(T):
			T = int(T)
		try:
			S_k = np.loadtxt(path+'sfac_T'+str(T)+'.txt')
		except:
			continue
		SkPeak = np.vstack([SkPeak, [T, S_k[k0]]])

	SkPeak[:,1] /= SkPeak[:,1].max()

	np.savetxt(path+'sfacpeak.txt', SkPeak, delimiter='\t', newline='\n')

	ax.plot(SkPeak[:,0], SkPeak[:,1], marker=markers[l], c=colors[l], markerfacecolor='None', lw=0.8, markersize=6, label=str(L))

ax.set_xlim(0,5)
ax.set_ylim(0,1.05)
ax.set_aspect(5/1.05)
ax.tick_params(direction='in', labelsize=22)
ax.set_xlabel('$T/J$', fontsize=24)
ax.set_ylabel('$S_{\\vec{q}\,=\,0}$', fontsize=24)
ax.set_xticks(np.arange(0,6,1))
ax.set_title('2D Ising FM', fontsize=24)
leg = ax.legend(title='$L$', fontsize=20, frameon=0)
leg.get_title().set_fontsize(24)
fig.tight_layout()
plt.show()
