import matplotlib.pyplot as plt
import numpy as np

def plot_simulation(u_prop, u_remi, u_nore, u_rocu, cp_prop, cp_remi, c_nore, cp_rocu, ce_prop, ce_remi, wav, bis, map,
                    co, hr, sv, tpr, nmb_m0, nmb_m1, nmb_m2, nmb_m3, titles: list = None, labels=None, disturbance_starts=None, volume_status_changes=None, use_cmap=False,
                    save_path=None):
    """
    Plot the simulation results.
    :param u_prop: Profofol infusion rate
    :param u_remi: Remifentanil infusion rate
    :param u_nore: Norepinephrine infusion rate
    :param u_rocu: Rocuronium infusion rate
    :param cp_prop: Propofol plasma concentration
    :param cp_remi: Remifentanil plasma concentration
    :param c_nore: Norepinephrine plasma concentration
    :param cp_rocu: Rocuronium plasma concentration
    :param ce_prop: Propofol effect site concentration
    :param ce_remi: Remifentanil effect site concentration
    :param wav: WAV index
    :param bis: BIS index
    :param map: Mean arterial pressure
    :param co: Cardiac output
    :param hr: Heart rate
    :param sv: Stroke volume
    :param nmb_m0: Neuromuscular blockade probability for m=0
    :param nmb_m1: Neuromuscular blockade probability for m=1
    :param nmb_m2: Neuromuscular blockade probability for m=2
    :param nmb_m3: Neuromuscular blockade probability for m=3
    :param titles: List of titles for the plots
    :param labels: List of labels for the plots
    :param use_cmap: Boolean to use colormap
    :param save_path: Path to save the plots
    :return: Figures of the simulation results [DoH, Hemodynamics, NMB]
    :rtype: tuple(matlplotlib.figure.Figure)
    """

    # Check if the input is a list of lists otherwise convert it to a list of lists to iterate over it
    if not all(isinstance(i, list) for i in u_prop):
        t_sim = [np.arange(0, len(u_prop))/60] # Convert to minutes
        u_prop, u_remi, u_nore, u_rocu, cp_prop, cp_remi, c_nore, cp_rocu, ce_prop, ce_remi, wav, bis, map, co, hr, sv,tpr, nmb_m0, nmb_m1, nmb_m2, nmb_m3 = [
            [x] for x in
            [u_prop, u_remi, u_nore, u_rocu, cp_prop, cp_remi, c_nore, cp_rocu, ce_prop, ce_remi, wav, bis, map, co, hr,
             sv, tpr, nmb_m0, nmb_m1, nmb_m2, nmb_m3]]
    else:
        t_sim = [np.arange(0, len(u))/60  for u in u_prop]

    if disturbance_starts:
        disturbance_starts = [disturbance_start/60 for disturbance_start in disturbance_starts]
    if volume_status_changes:
        volume_status_changes = [volume_status_change/60 for volume_status_change in volume_status_changes]

    # Helper function to get min and max from a list of lists
    def get_min_max(data):
        y_min  = min(min(values) for values in data) if isinstance(data, list) else min(data)
        y_max = max(max(values) for values in data) if isinstance(data, list) else max(data)
        return y_min * 0.8, y_max * 1.1

    # Calculate minimum and maximum values for each dataset
    y_min_bis, y_max_bis = get_min_max(bis)
    y_min_wav, y_max_wav = get_min_max(wav)
    y_min_map, y_max_map = get_min_max(map)
    y_min_co, y_max_co = get_min_max(co)
    y_min_hr, y_max_hr = get_min_max(hr)
    y_min_sv, y_max_sv = get_min_max(sv)
    y_min_tpr, y_max_tpr = get_min_max(tpr)

    # Create a colormap and a normalization instance if use_cmap is True
    cmap, norm = (plt.cm.viridis_r, plt.Normalize(vmin=0, vmax=len(u_prop) - 1)) if use_cmap else (None, None)

    def add_vertical_lines(ax, x_values, y_min, y_max, i, color='r', label=None):
        """Adds vertical lines to the axes"""
        if x_values and i == 0:
            # Only assign label to the first line to avoid duplicate legend entries
            for i, x in enumerate(x_values):
                current_label = label if i == 0 else None
                ax.axvline(x=x, color=color, linestyle='--', label=current_label),
        ax.set_ylim(y_min, y_max)

    def plot_data(ax, x, y,color, label, xlabel, ylabel, hlines_stimuli=None, hlines_volume_status=None, i=None,y_min=None,y_max=None):
        """Plots the data and adds supplementary vertical lines."""
        linestyle = '--' if 'Reference' in label else '-'
        ax.plot(x, y, color=color, label=label, linestyle=linestyle)
        # Add vertical lines
        if hlines_stimuli:
            add_vertical_lines(ax, hlines_stimuli, y_min, y_max,i ,color='r', label='Disturbance')
        if hlines_volume_status:
            add_vertical_lines(ax, hlines_volume_status, y_min, y_max,i, color='g', label='Volume Status Change')
        ax.set_xlabel(xlabel, fontsize=18)
        ax.set_ylabel(ylabel, fontsize=18)
        # Only draw legend if there are labeled lines
        if label or hlines_stimuli or hlines_volume_status:
            ax.legend(fontsize='xx-large')
        ax.grid(True)


    fig1, axs = plt.subplots(4, 2, figsize=(14, 16))
    fig1.suptitle(titles[0] if titles else f'DoH Simulation', fontsize=18)

    for i, (u_prop_val,u_remi_val,cp_prop_val, cp_remi_val, ce_prop_val, ce_remi_val, wav_val, bis_val) in enumerate(
            zip(u_prop,u_remi,cp_prop, cp_remi, ce_prop, ce_remi, wav, bis)):
        color = cmap(norm(i)) if cmap else None
        label = labels[i] if labels else ""
        plot_data(axs[0, 0], t_sim[i], u_prop_val, color, label, 'Time [min]', '$u_{Prop} [mg/s]$')
        plot_data(axs[0, 1], t_sim[i], u_remi_val, color, label, 'Time [min]', '$u_{Remi} [\mu g/s]$')
        plot_data(axs[1, 0], t_sim[i], cp_prop_val, color, label, 'Time [min]', '$C_{p, prop} [\mu g/ml]$')
        plot_data(axs[1, 1], t_sim[i], cp_remi_val, color, label, 'Time [min]', '$C_{p,remi} [ng/ml]$')
        plot_data(axs[2, 0], t_sim[i], ce_prop_val, color, label, 'Time [min]', '$C_{e,prop} [\mu g/ml]$')
        plot_data(axs[2, 1], t_sim[i], ce_remi_val, color, label, 'Time [min]', '$C_{e,remi} [ng/ml]$')
        plot_data(axs[3, 0], t_sim[i], wav_val, color, label, 'Time [min]', 'WAV$_{cns}$', disturbance_starts,volume_status_changes,i,y_min_wav,y_max_wav)
        plot_data(axs[3, 1], t_sim[i], bis_val, color, label, 'Time [min]', 'BIS', disturbance_starts,volume_status_changes,i,y_min_bis,y_max_bis)

    fig1.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig1.show()
    if save_path:
        title = titles[0] if titles else 'DoH Simulation'
        fig1.savefig(save_path / f"{title}.pdf", format='pdf', bbox_inches='tight')

    fig2, axs = plt.subplots(6, 2, figsize=(14, 24))
    fig2.suptitle(titles[1] if titles else 'Hemodynamic Variables Simulation', fontsize=16)

    for i, (u_prop_val, u_remi_val, u_nore_val, cp_prop_val, cp_remi_val, c_nore_val, map_val, co_val, hr_val,
            sv_val, tpr_val) in enumerate(zip(u_prop, u_remi, u_nore, cp_prop, cp_remi, c_nore, map, co, hr, sv, tpr)):
        color = cmap(norm(i)) if cmap else None
        label = labels[i] if labels else ""
        plot_data(axs[0, 0], t_sim[i], u_prop_val, color, label, 'Time [min]', '$u_{Prop} [mg/s]$')
        plot_data(axs[0, 1], t_sim[i], u_remi_val, color, label, 'Time [min]', '$u_{Remi} [\mu g/s]$')
        plot_data(axs[2, 0], t_sim[i], u_nore_val, color, label, 'Time [min]', '$u_{Nore} [\mu g/s]$')
        plot_data(axs[1, 0], t_sim[i], cp_prop_val, color, label, 'Time [min]', '$C_{p, prop} [\mu g/ml]$')
        plot_data(axs[1, 1], t_sim[i], cp_remi_val, color, label, 'Time [min]', '$C_{p,remi} [ng/ml]$')
        plot_data(axs[2, 1], t_sim[i], c_nore_val, color, label, 'Time [min]', '$C_{p,nore} [ng/ml]$')
        plot_data(axs[3, 0], t_sim[i], map_val, color, label, 'Time [min]', 'MAP [mmHg]',disturbance_starts,volume_status_changes,i,y_min_map,y_max_map)
        plot_data(axs[3, 1], t_sim[i], co_val, color, label, 'Time [min]', 'CO [L/min]',disturbance_starts,volume_status_changes,i,y_min_co ,y_max_co)
        plot_data(axs[4, 0], t_sim[i], hr_val, color, label, 'Time [min]', 'HR [beats/min]',disturbance_starts,volume_status_changes,i,y_min_hr,y_max_hr)
        plot_data(axs[4, 1], t_sim[i], sv_val, color, label, 'Time [min]', 'SV [ml]',disturbance_starts,volume_status_changes,i,y_min_sv,y_max_sv)
        plot_data(axs[5, 0], t_sim[i], tpr_val, color, label, 'Time [min]', 'TPR',disturbance_starts,volume_status_changes,i,y_min_tpr,y_max_tpr)
    axs[5,1].set_visible(False)

    fig2.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig2.show()
    if save_path:
        title = titles[1] if titles else 'Hemodynamic Variables Simulation'
        fig2.savefig(save_path / f"{title}.pdf", format='pdf', bbox_inches='tight')

    fig3, axs = plt.subplots(6, 1, figsize=(7, 24))
    fig3.suptitle(titles[2] if titles else 'Neuromuscular Blockade Simulation', fontsize=16)

    for i, (u_rocu_val, cp_rocu_val, nmb_m0_val, nmb_m1_val, nmb_m2_val, nmb_m3_val) in enumerate(
            zip(u_rocu, cp_rocu, nmb_m0, nmb_m1, nmb_m2, nmb_m3)):
        color = cmap(norm(i)) if cmap else None
        label = labels[i] if labels else ""
        plot_data(axs[0], t_sim[i], u_rocu_val, color, label, 'Time [min]', '$u_{Rocu} [mg/s]$')
        plot_data(axs[1], t_sim[i], cp_rocu_val, color, label, 'Time [min]', '$C_{p,Ro} [\mu g/ml]$')
        plot_data(axs[2], t_sim[i], nmb_m0_val, color, label, 'Time [min]', 'Probability (NMB m0)')
        plot_data(axs[3], t_sim[i], nmb_m1_val, color, label, 'Time [min]', 'Probability (NMB m1)')
        plot_data(axs[4], t_sim[i], nmb_m2_val, color, label, 'Time [min]', 'Probability (NMB m2)')
        plot_data(axs[5], t_sim[i], nmb_m3_val, color, label, 'Time [min]', 'Probability (NMB m3)')

    fig3.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig3.show()
    if save_path:
        title = titles[2] if titles else 'Neuromuscular Blockade Simulation'
        fig3.savefig(save_path / f"{title}.pdf", format='pdf', bbox_inches='tight')

    plt.close(fig1)
    plt.close(fig2)
    plt.close(fig3)

    return fig1, fig2, fig3