import sys
import h5py
import numpy as np
import matplotlib.pyplot as plt


f = h5py.File(sys.argv[1],"r")

t = f["Populations"]["GC"]["Cell Clamp Results"]["MPP ('AMPA', 'NMDA') PSC t"]["Attribute Value"][:]
v = f["Populations"]["GC"]["Cell Clamp Results"]["MPP ('AMPA', 'NMDA') PSC v"]["Attribute Value"][:]
i = f["Populations"]["GC"]["Cell Clamp Results"]["MPP ('AMPA', 'NMDA') PSC i mean"]["Attribute Value"][:]
i_var = f["Populations"]["GC"]["Cell Clamp Results"]["MPP ('AMPA', 'NMDA') PSC i variance"]["Attribute Value"][:]
i_AMPA = f["Populations"]["GC"]["Cell Clamp Results"]["MPP ('AMPA', 'NMDA') PSC i_AMPA"]["Attribute Value"][:]
i_NMDA = f["Populations"]["GC"]["Cell Clamp Results"]["MPP ('AMPA', 'NMDA') PSC i_NMDA"]["Attribute Value"][:]

vclamp_t = f["Populations"]["GC"]["Cell Clamp Results"]["MPP ('AMPA', 'NMDA') vclamp t"]["Attribute Value"][:]
ik = f["Populations"]["GC"]["Cell Clamp Results"]["MPP ('AMPA', 'NMDA') vclamp ik mean"]["Attribute Value"][:]
ina = f["Populations"]["GC"]["Cell Clamp Results"]["MPP ('AMPA', 'NMDA') vclamp ina mean"]["Attribute Value"][:]
i_Kir = f["Populations"]["GC"]["Cell Clamp Results"]["MPP ('AMPA', 'NMDA') vclamp i_Kir21 mean"]["Attribute Value"][:]
f.close()


vclamp_start_indexes = np.argwhere(vclamp_t == np.min(vclamp_t))
vclamp_parts = np.split(np.arange(len(vclamp_t)), np.asarray(vclamp_start_indexes.flat)[1:])

start_indexes = np.argwhere(t == np.min(t))
parts = np.split(np.arange(len(t)), np.asarray(start_indexes.flat)[1:])

V_holds = []
for n, part in enumerate(parts):
    vclamp_part = vclamp_parts[n]
    iclamp_period_idxs = np.argwhere(np.logical_and(t[part] > 200.0, t[part] < 450.0)).flat
    V_hold = np.round(np.mean(v[part][iclamp_period_idxs]))
    V_holds.append(V_hold)

fig, axs = plt.subplots(nrows=len(parts), ncols=1, sharex='col')
for n, part in enumerate(parts):
    vclamp_part = vclamp_parts[n]
    axs[n].plot(t[part], i[part], label="i")
    axs[n].fill_between(t[part], i[part] - i_var[part], i[part] + i_var[part])
    axs[n].plot(vclamp_t[vclamp_part], ik[vclamp_part], label="ik")
    axs[n].plot(vclamp_t[vclamp_part], ina[vclamp_part], label="ina")
    axs[n].title.set_text(f"V_hold = {V_holds[n]}")
    
plt.legend()
plt.show()

fig, axs = plt.subplots(nrows=len(parts), ncols=1, sharex='col')
for n, part in enumerate(parts):
    vclamp_part = vclamp_parts[n]
    axs[n].plot(t[part], i[part], label="i")
    axs[n].fill_between(t[part], i[part] - i_var[part], i[part] + i_var[part])
    axs[n].plot(vclamp_t[vclamp_part], i_Kir[vclamp_part], label="i_Kir")
    axs[n].plot(vclamp_t[vclamp_part], ina[vclamp_part], label="ina")
    
plt.legend()
plt.show()


fig, axs = plt.subplots(nrows=len(parts), ncols=1, sharex='col')
for n, part in enumerate(parts):
    axs[n].plot(t[part], i_AMPA[part], label='AMPA')
    axs[n].plot(t[part], i_NMDA[part], label='NMDA')

plt.legend()
plt.show()
