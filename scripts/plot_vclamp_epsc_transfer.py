import sys
import h5py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from scipy import integrate

font = {'family' : 'sans-serif',
        'sans-serif': 'Arial',
        'style': 'normal',
        'weight': 'normal',
        'size'   : 24 }
matplotlib.rc('font', **font)


t1 = 200.0
t2 = 450.0

dfs = []

grp = "NGFC ('GABA_A', 'GABA_B')"
grp = "MC ('AMPA', 'NMDA')"
grp = "LPP ('AMPA', 'NMDA')"
grp = "MPP ('AMPA', 'NMDA')"

for fname in sys.argv[1:]:
    f = h5py.File(fname,"r")

    t = f["Populations"]["GC"]["Cell Clamp Results"][f"{grp} PSC t"]["Attribute Value"][:]
    v = f["Populations"]["GC"]["Cell Clamp Results"][f"{grp} PSC v"]["Attribute Value"][:]
    i = f["Populations"]["GC"]["Cell Clamp Results"][f"{grp} PSC i mean"]["Attribute Value"][:]
#    i_AMPA = f["Populations"]["GC"]["Cell Clamp Results"][f"{grp} PSC i_AMPA"]["Attribute Value"][:]
#    i_NMDA = f["Populations"]["GC"]["Cell Clamp Results"][f"{grp} PSC i_NMDA"]["Attribute Value"][:]
    
    vclamp_t = f["Populations"]["GC"]["Cell Clamp Results"][f"{grp} vclamp t"]["Attribute Value"][:]
    ik = f["Populations"]["GC"]["Cell Clamp Results"][f"{grp} vclamp ik mean"]["Attribute Value"][:]
    ik_var = f["Populations"]["GC"]["Cell Clamp Results"][f"{grp} vclamp ik variance"]["Attribute Value"][:]
    ina = f["Populations"]["GC"]["Cell Clamp Results"][f"{grp} vclamp ina mean"]["Attribute Value"][:]
    ina_var = f["Populations"]["GC"]["Cell Clamp Results"][f"{grp} vclamp ina variance"]["Attribute Value"][:]
    i_Kir = f["Populations"]["GC"]["Cell Clamp Results"][f"{grp} vclamp ik_Kir21 mean"]["Attribute Value"][:]
    i_Kir_var = f["Populations"]["GC"]["Cell Clamp Results"][f"{grp} vclamp ik_Kir21 variance"]["Attribute Value"][:]
    i_KA = f["Populations"]["GC"]["Cell Clamp Results"][f"{grp} vclamp ik_Kv42 mean"]["Attribute Value"][:]
    i_KA_var = f["Populations"]["GC"]["Cell Clamp Results"][f"{grp} vclamp ik_Kv42 variance"]["Attribute Value"][:]
    i_fKDR = f["Populations"]["GC"]["Cell Clamp Results"][f"{grp} vclamp ik_Kv11 mean"]["Attribute Value"][:]
    i_fKDR_var = f["Populations"]["GC"]["Cell Clamp Results"][f"{grp} vclamp ik_Kv11 variance"]["Attribute Value"][:]
    i_sKDR = None
    i_sKDR_var = None
    if f"{grp} vclamp i_sKDR_Aradi mean" in f["Populations"]["GC"]["Cell Clamp Results"]:
        i_sKDR = f["Populations"]["GC"]["Cell Clamp Results"][f"{grp} vclamp i_sKDR_Aradi mean"]["Attribute Value"][:]
        i_sKDR_var = f["Populations"]["GC"]["Cell Clamp Results"][f"{grp} vclamp i_sKDR_Aradi variance"]["Attribute Value"][:]
    f.close()

    vclamp_start_indexes = np.argwhere(vclamp_t == np.min(vclamp_t))
    vclamp_parts = np.split(np.arange(len(vclamp_t)), np.asarray(vclamp_start_indexes.flat)[1:])

    start_indexes = np.argwhere(t == np.min(t))
    parts = np.split(np.arange(len(t)), np.asarray(start_indexes.flat)[1:])

    V_holds = []
    for n, part in enumerate(parts):
        vclamp_part = vclamp_parts[n]
        vclamp_period_idxs = np.argwhere(np.logical_and(t[part] >= t1, t[part] <= t2)).flat
        V_hold = np.round(np.mean(v[part][vclamp_period_idxs]))
        V_holds.append(V_hold)
    
    i_ints = []
    ik_ints = []
    i_Kir_ints = []
    i_KA_ints = []
    i_fKDR_ints = []
    i_sKDR_ints = []
    ina_ints = []
    for n, part in enumerate(parts):
        vclamp_part = vclamp_parts[n]
        iclamp_period_idxs = np.argwhere(np.logical_and(t[part] >= t1, t[part] <= t2)).flat
        vclamp_period_idxs = np.argwhere(np.logical_and(vclamp_t[vclamp_part] >= t1, vclamp_t[vclamp_part] <= t2)).flat

        i_int = integrate.simpson(i[part][iclamp_period_idxs], x=t[part][iclamp_period_idxs])
        ik_int = integrate.simpson(ik[vclamp_part][vclamp_period_idxs], x=vclamp_t[vclamp_part][vclamp_period_idxs])
        ina_int = integrate.simpson(ina[vclamp_part][vclamp_period_idxs], x=vclamp_t[vclamp_part][vclamp_period_idxs])
        i_Kir_int = integrate.simpson(i_Kir[vclamp_part][vclamp_period_idxs], x=vclamp_t[vclamp_part][vclamp_period_idxs])
        i_KA_int = integrate.simpson(i_KA[vclamp_part][vclamp_period_idxs], x=vclamp_t[vclamp_part][vclamp_period_idxs])
        i_fKDR_int = integrate.simpson(i_fKDR[vclamp_part][vclamp_period_idxs], x=vclamp_t[vclamp_part][vclamp_period_idxs])
        if i_sKDR is not None:
            i_sKDR_int = integrate.simpson(i_sKDR[vclamp_part][vclamp_period_idxs], x=vclamp_t[vclamp_part][vclamp_period_idxs])
    
        i_ints.append(i_int)
        ik_ints.append(ik_int)
        ina_ints.append(ina_int)
        i_Kir_ints.append(i_Kir_int)
        i_KA_ints.append(i_KA_int)
        i_fKDR_ints.append(i_fKDR_int)
        if i_sKDR is not None:
            i_sKDR_ints.append(i_sKDR_int)

    q_dict = { 'i_int': i_ints,
               'ik_int': ik_ints,
               'ina_int': ina_ints,
               'i_Kir_int': i_Kir_ints,
               'i_KA_int': i_KA_ints,
               'i_fKDR_int': i_fKDR_ints,
               'V_hold': V_holds }
    if i_sKDR is not None:
        q_dict['i_sKDR_int'] = i_sKDR_ints
    df = pd.DataFrame.from_dict(q_dict)
    dfs.append(df)


df = pd.concat(dfs)

V_hold_levels = df["V_hold"].unique()

int_amp_mean_vars = []

q_names = [ "i_int", "ik_int", "ina_int", "i_int_amp",
            "i_Kir_int", "i_KA_int", "i_fKDR_int", ]
if i_sKDR is not None:
    q_names.append("i_sKDR_int")
for V_hold in V_hold_levels:

    sub_df = df.loc[df["V_hold"] == V_hold]
    i_int_amp = sub_df[["i_int", "ik_int", "ina_int"]].sum(axis=1)
    sub_df.insert(loc=0, column='i_int_amp', value=i_int_amp)
    q_means = []
    q_vars  = []
    
    for q in q_names:
        mean = sub_df[q].mean()
        var = sub_df[q].var()
        q_means.append(mean)
        q_vars.append(var)
    
    int_amp_mean_vars.append(pd.DataFrame.from_dict( { 'quantity': q_names,
                                                        'V_hold': np.asarray([V_hold]*len(q_names)),
                                                        'mean': q_means,
                                                        'variance': q_vars }))
int_amp_mean_var_df = pd.concat(int_amp_mean_vars)

fig, axs = plt.subplots(nrows=1, ncols=2)

sub_df = int_amp_mean_var_df.loc[int_amp_mean_var_df["quantity"] == "i_int_amp"]
y = sub_df["mean"]
error = sub_df["variance"]
axs[0].plot(np.asarray(V_holds), y, linewidth=3)
axs[0].fill_between(np.asarray(V_holds), y-error, y+error)
axs[0].set_ylabel("Charge Transfer [uC]")
axs[0].set_xlabel("Voltage [mV]")

display_q_names = ["i_int", "ik_int", "ina_int",
                   "i_Kir_int", "i_KA_int", "i_fKDR_int"]

if i_sKDR is not None:
    display_q_names.append("i_sKDR_int")

for q in display_q_names:
    sub_df = int_amp_mean_var_df.loc[int_amp_mean_var_df["quantity"] == q]
    y = sub_df["mean"]
    error = sub_df["variance"]
    axs[1].plot(np.asarray(V_holds), y, label=q, linewidth=3)
    axs[1].fill_between(np.asarray(V_holds), y-error, y+error)
    
axs[1].set_xlabel("Voltage [mV]")
    
plt.legend(fontsize='x-small')
plt.show()
    
