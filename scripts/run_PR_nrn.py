import os, sys, logging, click, yaml, pprint
import numpy as np
from neuron import h
from scipy import signal

import matplotlib.pyplot as plt

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def list_find(f, lst):
    """

    :param f:
    :param lst:
    :return:
    """
    i = 0
    for x in lst:
        if f(x):
            return i
        else:
            i = i + 1
    return None


def detect_spikes(T, Y, t0, t1, before_peak=50.0):

    spk_info_dtype = np.dtype(
        [
            ("Vpeak", float),
            ("Tpeak", float),
            ("amplitude", float),
            ("T0", float),
            ("T1", float),
        ]
    )

    dt = np.mean(np.diff(T))
    pre_period_idxs = np.argwhere(T < t0 - before_peak).flat
    pre_peak_info = signal.find_peaks(
        Y[pre_period_idxs], height=-20.0, width=(None, int(before_peak / dt))
    )
    pre_peak_idxs = pre_peak_info[0]
    N_peaks_pre = len(pre_peak_idxs)

    spk_period_idxs = np.argwhere(np.logical_and(T >= t0 - before_peak, T <= t1)).flat
    T_spk = T[spk_period_idxs]
    Y_spk = Y[spk_period_idxs]
    peak_info = signal.find_peaks(
        Y_spk, height=-20.0, width=(None, int(before_peak / dt))
    )
    peak_idxs = peak_info[0]
    if len(peak_idxs) == 0:
        return N_peaks_pre, 0, None, None

    # Determine threshold based on the following:
    # 1. take the dV/dt of the voltage trace
    # 2. measure the SD of the dV/dt for 50 ms before the AP.
    # 3. threshold = mean(dV/dt) + 2*SD(dV/dt) (Atherton and Bevan 2005).
    dydt = np.gradient(Y_spk, T_spk)
    peak_idx = peak_idxs[0]
    T_peak = T_spk[peak_idx]
    T_before_idxs = np.argwhere(
        np.isclose(T_spk, T_peak - before_peak, rtol=1e-4, atol=1e-4)
    )
    if len(T_before_idxs) == 0:
        return N_peaks_pre, 0, None, None

    T_before_idx = T_before_idxs[0][0]
    mean_dydt = np.mean(dydt[T_before_idx:peak_idx])
    sd_dydt = np.std(dydt[T_before_idx:peak_idx])
    dydt_threshold = mean_dydt + 2 * sd_dydt
    try:
        threshold_idx = np.argmin(np.abs(dydt[T_before_idx:peak_idx] - dydt_threshold))
    except:
        logger.debug(
            f"Error in threshold computation: dydt_threshold = {dydt_threshold} dydt[T_before_idx:peak_idx] = {dydt[T_before_idx:peak_idx]} T[T_before_idx:peak_idx] = {T[T_before_idx:peak_idx]} Y[T_before_idx:peak_idx] = {Y[T_before_idx:peak_idx]} "
        )
        return N_peaks_pre, 0, None, None

    threshold = Y_spk[T_before_idx:peak_idx][threshold_idx]

    period_idxs = np.argwhere(np.logical_and(T >= t0, T <= t1)).flat
    T = T[period_idxs]
    Y = Y[period_idxs]

    # 1. Make the data binary, in a way that they are true when larger than the threshold and false when lower or equal.
    # 2. Take the difference of the binary signal.
    threshold_crossings = np.diff(Y > threshold, prepend=False)
    crossing_idx = np.argwhere(threshold_crossings)[:, 0]
    up_crossing_idx = np.argwhere(threshold_crossings)[::2, 0]
    N_peaks = len(up_crossing_idx)

    spk_info = None
    peak_amps = None
    T_peaks = None
    Y_peaks = None
    if N_peaks > 0:

        # Split V and T into intervals based on threshold crossing indices
        Y_intervals = np.split(Y, crossing_idx[1::2])[:-1]
        T_intervals = np.split(T, crossing_idx[1::2])[:-1]
        # Obtain peak indices in each V interval
        peak_idxs = [np.argmax(Y_interval) for Y_interval in Y_intervals]
        N_peaks = len(peak_idxs)
        Y_peaks = []
        T_peaks = []
        peak_amps = []
        for j, (T_interval, peak_idx) in enumerate(zip(T_intervals, peak_idxs)):
            if len(T_interval) < 2:
                N_peaks -= 1
                continue
            peak_amp = np.max(Y_intervals[j]) - threshold
            peak_amps.append(peak_amp)
            Y_peak = Y_intervals[j][peak_idx]
            Y_peaks.append(Y_peak)
            if peak_idx < len(T_interval):
                T_peaks.append(T_interval[peak_idx])
            else:
                T_peaks.append(T_interval[-1])

        if N_peaks > 0:
            spk_info = np.zeros(shape=(N_peaks), dtype=spk_info_dtype)
            for p in range(N_peaks):
                spk_info[p]["Vpeak"] = Y_peaks[p]
                spk_info[p]["Tpeak"] = T_peaks[p]
                spk_info[p]["T0"] = T_intervals[p][0]
                spk_info[p]["T1"] = T_intervals[p][-1]
                spk_info[p]["amplitude"] = peak_amps[p]

    return N_peaks_pre, N_peaks, threshold, spk_info

def run_iclamp(
    cell,
    amp,
    t0,
    t1,
    dt=0.025,
    record_dt=0.01,
    t_stop=1000.0,
    v_init=-65.0,
    celsius=36,
):

    # Create the recording vectors for time and voltage
    vec_t = h.Vector()
    vec_soma_v = h.Vector()
    vec_dend_v = h.Vector()
    vec_dend_ica = h.Vector()
    vec_dend_ik = h.Vector()
    vec_soma_ik = h.Vector()
    vec_soma_ina = h.Vector()
    vec_dend_cai = h.Vector()
    vec_dend_ki = h.Vector()
    vec_soma_ki = h.Vector()
    vec_soma_g_Na = h.Vector()
    vec_dend_g_KAHP = None
    vec_dend_g_KCa = h.Vector()
    vec_dend_g_Ca = h.Vector()
    vec_dend_g_HCN = None

    vec_t.record(h._ref_t, record_dt)  # Time
    vec_soma_v.record(cell.soma(0.5)._ref_v, record_dt)  # Voltage
    vec_dend_v.record(cell.dend(0.5)._ref_v, record_dt)  # Voltage
    vec_dend_ica.record(cell.dend(0.5)._ref_ica, record_dt)
    vec_dend_ik.record(cell.dend(0.5)._ref_ik, record_dt)
    vec_soma_ik.record(cell.soma(0.5)._ref_ik, record_dt)
    vec_soma_ina.record(cell.soma(0.5)._ref_ina, record_dt)
    vec_dend_cai.record(cell.dend(0.5)._ref_cai, record_dt)
    vec_dend_ki.record(cell.dend(0.5)._ref_ki, record_dt)
    vec_soma_ki.record(cell.soma(0.5)._ref_ki, record_dt)

    vec_soma_g_Na.record(cell.soma(0.5)._ref_g_Na_PR, record_dt)
    vec_dend_g_KCa.record(cell.dend(0.5)._ref_g_KCa_PR, record_dt)
    if h.ismembrane("KAHP_PR", sec=cell.dend):
        vec_dend_g_KAHP = h.Vector()
        vec_dend_g_KAHP.record(cell.dend(0.5)._ref_g_KAHP_PR, record_dt)
    vec_dend_g_Ca.record(cell.dend(0.5)._ref_g_Ca_PR, record_dt)
    if h.ismembrane("HCN", sec=cell.dend):
        vec_dend_g_HCN = h.Vector()
        vec_dend_g_HCN.record(cell.dend(0.5)._ref_g_HCN, record_dt)

    # Put an IClamp at the soma
    stim = h.IClamp(0.5, sec=cell.soma)
    stim.delay = t0  # Stimulus stat
    stim.dur = t1 - t0  # Stimulus length
    stim.amp = amp  # strength of current injection

    # Run the Simulation
    h.dt = dt
    h.celsius = celsius
    h.v_init = v_init
    h.init()
    h.finitialize(h.v_init)
    logger.info(f"soma ek = {cell.soma.ek} ena = {cell.soma.ena}")
    logger.info(f"dend ek = {cell.dend.ek} eca = {cell.dend.eca}")
    logger.info(f"soma nao = {cell.soma.nao} ko = {cell.soma.ko}")
    logger.info(f"dend ki = {cell.soma.ki} ko = {cell.dend.ko}")
    logger.info(
        f"dend cao = {cell.dend.cao} cai = {cell.dend.cai}"
    )

    h.tstop = t_stop
    h.run()

    result_dict = {
        "t": np.array(vec_t),
        "soma_v": np.array(vec_soma_v),
        "dend_v": np.array(vec_dend_v),
        "soma_ik": np.array(vec_soma_ik),
        "soma_ki": np.array(vec_soma_ki),
        "soma_ina": np.array(vec_soma_ina),
        "dend_ica": np.array(vec_dend_ica),
        "dend_cai": np.array(vec_dend_cai),
        "dend_ki": np.array(vec_dend_ki),
        "dend_ik": np.array(vec_dend_ik),
        "soma_g_Na": np.array(vec_soma_g_Na),
        "dend_g_Ca": np.array(vec_dend_g_Ca),
        "dend_g_KCa": np.array(vec_dend_g_KCa),
    }
    if vec_dend_g_KAHP is not None:
        result_dict["dend_g_KAHP"] = vec_dend_g_KAHP
    if vec_dend_g_HCN is not None:
        result_dict["dend_g_HCN"] = vec_dend_g_HCN

    return result_dict


@click.command()
@click.option(
    "--config-path",
    "-c",
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    help="path to configuration file",
)
@click.option("--model-variant", "-m", default="default", type=str)
@click.option("--dt", type=float, default=0.025, help="default simulation time step")
@click.option("--cvode/--no-cvode", default=False, help="use adaptive time step solver")
@click.option(
    "--stim-amp", type=float, default=0.08, help="amplitude of stimulus current"
)
@click.option("--stim-start", type=float, default=500.0, help="start time of stimulus")
@click.option("--stim-stop", type=float, default=1000.0, help="stop time of stimulus")
@click.option(
    "--t-stop", "-t", type=float, default=2000.0, help="stop time of simulation"
)
@click.option(
    "--template-path",
    required=False,
    default="templates",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    help="path to directory containing template file",
)
@click.option(
    "--template-file",
    required=False,
    default="PR_nrn.hoc",
    type=str,
    help="name of template file",
)
@click.option(
    "--template-name",
    required=False,
    default="PR_nrn",
    type=str,
    help="name of template class",
)
@click.option(
    "--nrn-type",
    required=False,
    default="PinskyRinzel",
    type=str,
    help="name of neuron type",
)
@click.option(
    "--param-key",
    type=str,
    default="PinskyRinzel",
    help="key for toplevel section with parameters",
)
def main(
    config_path,
    model_variant,
    dt,
    cvode,
    stim_amp,
    stim_start,
    stim_stop,
    t_stop,
    template_path,
    template_file,
    template_name,
    nrn_type,
    param_key,
):

    # Load the NEURON libraries
    h.load_file("stdrun.hoc")
    h.load_file("rn.hoc")

    # Enable variable time step solver
    h.cvode.use_fast_imem(1)
    h.cvode.cache_efficient(1)
    h.cvode.active(1 if cvode else 0)
    h.secondorder = 2
    h.dt = dt

    config_dict = None
    with open(config_path) as f:
        config_dict = yaml.load(f, Loader=yaml.FullLoader)

    h.load_file(os.path.join(template_path, template_file))
    param_dict = config_dict.get(param_key, None)
    if param_dict is None:
        param_dict = config_dict[int(param_key)]
    logger.info(f"{pprint.pformat(param_dict)}")

    v_init = param_dict["V_rest"]

    template = getattr(h, template_name)
    cell = template(param_dict)

    h.v_init = v_init
    h.init()
    h.finitialize(h.v_init)
    cell.init_ic(h.v_init)
    cell.soma.ic_constant = param_dict["ic_constant"]

    h.finitialize(h.v_init)
    h.finitialize(h.v_init)

    h.psection(sec=cell.soma)
    h.psection(sec=cell.dend)


    iclamp_results = run_iclamp(
        cell,
        stim_amp,
        stim_start,
        stim_stop,
        t_stop=t_stop,
        v_init=v_init,
        dt=dt,
    )
    vec_t = iclamp_results["t"]
    vec_soma_v = iclamp_results["soma_v"]
    vec_dend_v = iclamp_results["dend_v"]
    vec_dend_ik = iclamp_results["dend_ik"]
    vec_soma_ik = iclamp_results["soma_ik"]
    vec_soma_ina = iclamp_results["soma_ina"]
    vec_dend_ik = iclamp_results["dend_ik"]
    vec_dend_ica = iclamp_results["dend_ica"]
    vec_dend_cai = iclamp_results["dend_cai"]
    vec_dend_ki = iclamp_results["dend_ki"]
    vec_soma_ki = iclamp_results["soma_ki"]

    vec_soma_g_Na = iclamp_results["soma_g_Na"]
    vec_dend_g_KAHP = iclamp_results.get("dend_g_KAHP", None)
    vec_dend_g_KCa = iclamp_results["dend_g_KCa"]
    vec_dend_g_Ca = iclamp_results["dend_g_Ca"]
    vec_dend_g_HCN = iclamp_results.get("dend_g_HCN", None)

    print(f"ic_constant: {cell.soma.ic_constant}")
    logger.info(f"spikes: {detect_spikes(vec_t, vec_soma_v, stim_start, stim_stop+5.0)}")
    logger.info(f"dend ki min/max: {np.min(vec_dend_ki)} / {np.max(vec_dend_ki)}")

    nrows = 6
    ncols = 3
    fig, axs = plt.subplots(nrows, ncols)
    axs[0, 0].plot(vec_t, vec_soma_v, linewidth=3, color="r", label="soma_v")
    axs[1, 0].plot(vec_t, vec_dend_v, linewidth=3, color="r", label="dend_v")
    axs[2, 0].plot(vec_t, vec_soma_ina, linewidth=3, color="b", label="soma_ina")
    axs[3, 0].plot(vec_t, vec_dend_ik, linewidth=3, color="b", label="soma_ik")
    axs[4, 0].plot(vec_t, vec_dend_ik, linewidth=3, color="b", label="dend_ik")
    axs[5, 0].plot(vec_t, vec_dend_ica, linewidth=3, color="r", label="dend_ica")
    axs[-1, 0].set_xlabel("Time (ms)")
    axs[0, 0].set_ylabel("V (mV)")

    axs[0, 1].plot(vec_t, vec_dend_g_Ca, linewidth=3, color="b", label="dend_g_Ca")
    axs[1, 1].plot(vec_t, vec_dend_g_KCa, linewidth=3, color="b", label="dend_g_KCa")
    if vec_dend_g_KAHP is not None:
        axs[2, 1].plot(vec_t, vec_dend_g_KAHP, linewidth=3, color="b", label="dend_g_KAHP")
    if vec_dend_g_HCN is not None:
        axs[3, 1].plot(vec_t, vec_dend_g_HCN, linewidth=3, color="b", label="dend_g_HCN")
        
    axs[0, 2].plot(vec_t, vec_dend_cai, linewidth=3, color="g", label="dend_cai")
    axs[1, 2].plot(vec_t, vec_dend_ki, linewidth=3, color="g", label="dend_ki")
    axs[2, 2].plot(vec_t, vec_soma_ki, linewidth=3, color="g", label="soma_ki")

    for i in range(nrows):
        for j in range(ncols):
            axs[i, j].legend()
    plt.savefig(f"{nrn_type}_iclamp.svg")
    plt.show()


if __name__ == "__main__":
    main(
        args=sys.argv[
            (
                list_find(
                    lambda x: os.path.basename(x) == os.path.basename(__file__),
                    sys.argv,
                )
                + 1
            ) :
        ]
    )
