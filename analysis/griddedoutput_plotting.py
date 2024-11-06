import matplotlib.pyplot as plt
from cycler import cycler
#from griddedoutput_loaddatasets import *
from loaddatastructs import * # use GriddedOutput data object
from griddedoutput_helperfuncs import *
plt.rcParams['axes.prop_cycle'] = cycler('color', plt.get_cmap('tab20').colors)

def _timeLabel(t_idx, ax, **kwargs):
    time = int((t_idx - 1)*GriddedOutput.historydelta_m/60)

    fontsize=kwargs.get('subplot_title_fontsize', 12)
    
    if time == 0:
        time_text = 'Initial condition'
    else:
        time_text = f'{time} hours'

    txt = AnchoredText(f'{time_text}',
                    loc='lower left', prop=dict(size=fontsize), frameon=True,
                    bbox_to_anchor=(0., 1.), borderpad=0, pad=0.25,
                    bbox_transform=ax.transAxes
                    )
    ax.add_artist(txt)

def plotNumberDist(number_dist_arr, plot_type='hist', **kwargs):
    if not kwargs.get('ax'):
        fig, ax = plt.subplots(1, 1, figsize=(8,6))
    else:
        ax = kwargs.get('ax')

    if plot_type == 'line':
        c = kwargs.get('color', 'k')
        label = kwargs.get('label', None)

        if kwargs.get('t_idx'):
            t = kwargs.get('t_idx')
            history_dt = GriddedOutput.historydelta_m/60 # hours
            time =  ((t-1)*history_dt)
            label = f't = {time:3.2f} hrs'

        lw = kwargs.get('lw', 1)
        ls = kwargs.get('ls', '-')
        ax.plot(GriddedOutput.bin_geocenter, number_dist_arr, c=c, label=label, lw=lw, ls=ls)
    elif plot_type == 'hist':
        bin_width = GriddedOutput.bin_edges[1:] - GriddedOutput.bin_edges[:-1] 
        ax.bar(x=GriddedOutput.bin_geocenter, height=number_dist_arr, width=bin_width, alpha=.6)
    if kwargs.get('ylims'):
        ax.set_ylim(kwargs.get('ylims'))
    if kwargs.get('xlims'):
        ax.set_xlim(kwargs.get('xlims'))
    ax.set_xscale('log')
    ax.set_yscale(kwargs.get('yscale', 'log'))
    ax.grid(which = "major", linewidth = 1)
    ax.grid(which = "minor", linewidth = 0.2)
    ax.minorticks_on()
    ax.set_ylabel('$dN/d\log{D_p}$')
    ax.set_xlabel('Particle Diameter [m]')
    leg_fontsize = kwargs.get('legend_fontsize', 10)

    ax.legend(fontsize=leg_fontsize)
    if kwargs.get('ax'):
        return ax

def plotSpeciatedMassDist(binned_species_mass_arr, **kwargs):

    if binned_species_mass_arr.shape[1] == GriddedOutput.n_bins:
        binned_species_mass_arr = binned_species_mass_arr[:, :-1]

    bin_total_mass = np.zeros(((GriddedOutput.n_bins-1),))
    fig, ax = plt.subplots(1, 1)
    plt.rcParams['axes.prop_cycle'] = cycler('color', plt.get_cmap('tab20').colors)

    # Plot the total speciated mass distribution
    for i, species in enumerate(GriddedOutput.aero_species):
        species_mass_dist = binned_species_mass_arr[i, :]#/bin_logwidth # moved division by logwidth to calculation of mass distrib
        bin_total_mass += species_mass_dist
        #print(bin_total_mass[50])
    #ax.plot(bin_edges[:-1], bin_total_mass, label=species)
        if i == 0:
            y_lower = ax.get_ylim()[0]
            y_lower = np.array((GriddedOutput.n_bins-1)*[y_lower])
        else:
            y_lower = bin_total_mass - species_mass_dist
        ax.fill_between(x=GriddedOutput.bin_geocenter, y1=y_lower, y2=bin_total_mass, label=species)
    
    plt.xscale('log')
    plt.yscale(kwargs.get('yscale', 'log'))
    plt.ylabel('Mass [kg$\cdot$m$^{-3}$]')
    plt.xlabel('Particle Diameter [m]')
    if kwargs.get('xlims'):
        plt.xlim(kwargs.get('xlims'))
    plt.legend(loc='center', bbox_to_anchor=(1.3, .5), ncol=2)
    if kwargs.get('t_idx'):
        t = kwargs.get('t_idx')
        history_dt = GriddedOutput.historydelta_m/60 # hours
        time =  ((t-1)*history_dt)
        plt.text(1.3, .9, f't = {time:3.2f} hrs', transform=ax.transAxes, horizontalalignment='center')

    if kwargs.get('savefig'):
        filename = construct_figure_filename(figure_type='speciated-mass-dist', **kwargs)
        plt.savefig(filename, format='pdf', bbox_inches='tight')

def plotSpeciatedMassFrac(binned_species_mass_arr, **kwargs):
    if binned_species_mass_arr.shape[1] == GriddedOutput.n_bins:
        binned_species_mass_arr = binned_species_mass_arr[:, :-1]

    total_mass_per_bin = (binned_species_mass_arr[:, :]#/bin_logwidth # moved division by logwidth to calculation of mass distrib
                          ).sum(axis=0)
    frac_total = np.zeros((GriddedOutput.n_bins-1,))
    fig, ax = plt.subplots(1, 1)
    # Plot the total speciated mass distribution
    for i, species in enumerate(GriddedOutput.aero_species):
        species_mass_dist = binned_species_mass_arr[i, :]#/bin_logwidth # moved division by logwidth to calculation of mass distrib
        species_frac = species_mass_dist/total_mass_per_bin
        frac_total += species_frac
        if i == 0:
            y_lower = 0
            y_lower = np.array((GriddedOutput.n_bins-1)*[y_lower])
        else:
            y_lower = frac_total - species_frac
        ax.fill_between(x=GriddedOutput.bin_geocenter, y1=y_lower, y2=frac_total, label=species)

    plt.xlim(GriddedOutput.bin_edges[0], GriddedOutput.bin_edges[-1])
    plt.ylim(0, 1)
    plt.xscale('log')
    plt.ylabel('Mass Fraction')
    plt.xlabel('Particle Diameter [m]')
    if kwargs.get('xlims'):
        plt.xlim(kwargs.get('xlims'))
    plt.legend(loc='center', bbox_to_anchor=(1.3, .5), ncol=2)
    if kwargs.get('t_idx'):
        t = kwargs.get('t_idx')
        history_dt = GriddedOutput.historydelta_m/60 # hours
        time =  ((t-1)*history_dt)
        plt.text(1.3, .9, f't = {time:3.2f} hrs', transform=ax.transAxes, horizontalalignment='center')

    if kwargs.get('savefig'):
        filename = construct_figure_filename(figure_type='speciated-mass-frac', **kwargs)
        plt.savefig(filename, format='pdf', bbox_inches='tight')

def construct_figure_filename(figure_type, **kwargs):
    if kwargs.get('scenario'):
        filename = kwargs.get('scenario')
        filename += f'_{figure_type}'
    else:
        filename = figure_type

    if kwargs.get('t_idx'):
        filename +=  f'_t{int(kwargs.get("t_idx"))}'
    if kwargs.get('z_idx'):
        filename +=  f'_z{int(kwargs.get("z_idx"))}'
    
    extension = kwargs.get('figure_extension', 'pdf')
    filename += f'.{extension}'
    return filename

def plotFourPanelMassFrac(scenario, times, xstart, xend, ystart, yend, z_idx, **kwargs):
    
    fig, axs  = plt.subplots(2,2, figsize=(6.5,4.5))
    plt.subplots_adjust(hspace=.3, wspace=.2)

    if len(times) != 4:
        raise AttributeError("Number of times must be four")
        
    for j, (ax, t) in enumerate(zip(axs.flatten(), times)):

        GriddedOutput.loadData(scenario, xstart, xend, ystart, yend, z_idx, t)
        binned_species_mass_arr =  getBinnedSpeciesMassOptimized(GriddedOutput.gridded_data[scenario]['aero_diams'],
                                                                GriddedOutput.gridded_data[scenario]['aero_masses'], 
                                                                GriddedOutput.gridded_data[scenario]['aero_numconc'], 
                                                                n_grid_cells=GriddedOutput.gridded_data[scenario]['n_total_cells'])


        if binned_species_mass_arr.shape[1] == GriddedOutput.n_bins:
            binned_species_mass_arr = binned_species_mass_arr[:, :-1]

        total_mass_per_bin = (binned_species_mass_arr[:, :]#/bin_logwidth # moved division by logwidth to calculation of mass distrib
                            ).sum(axis=0)
        frac_total = np.zeros((GriddedOutput.n_bins-1,))

        # Plot the total speciated mass distribution
        aero_species = [var for var in GriddedOutput.aero_vars if 'pmc_' in var]
        for i, species in enumerate(aero_species):
            species_mass_dist = binned_species_mass_arr[i, :]#/bin_logwidth # moved division by logwidth to calculation of mass distrib
            species_frac = species_mass_dist/total_mass_per_bin
            frac_total += species_frac
            if i == 0:
                y_lower = 0
                y_lower = np.array((GriddedOutput.n_bins-1)*[y_lower])
            else:
                y_lower = frac_total - species_frac

            variable_fmt = GriddedOutput.aerosol_fmt_map[species]
            variable_fmt = variable_fmt.replace('Aerosol ', '')
            if ('ARO' in species) or ('ALK' in species) or ('OLE' in species) or ('API' in species) or ('LIM' in species):
                variable_fmt = variable_fmt.replace('_', '')
            if variable_fmt in ['CL', 'NA', 'CA']:
                variable_fmt = variable_fmt.title()
            ax.fill_between(x=GriddedOutput.bin_geocenter, y1=y_lower, y2=frac_total, label=variable_fmt)

        # set xlim, ylim, axis scaling
        ax.set_xlim(GriddedOutput.bin_edges[0], GriddedOutput.bin_edges[-1])
        ax.set_ylim(0, 1)
        ax.set_xscale('log')
        if kwargs.get('xlims'):
            ax.set_xlim(kwargs.get('xlims'))
        
        # add legend
        if j == 0:
            fig.legend(fontsize=11, ncol=5, loc='center', bbox_to_anchor=(.5,-.15))

        # Set x-axis ticks and labels
        if j < 2:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel('Particle Diameter [m]', fontsize=11)

        # Set y-axis ticks and labels
        if j in [1,3]:
            ax.set_yticklabels([])
        else:
            # Set y-axis ticks and label
            ax.set_ylabel('Mass Fraction', fontsize=11)  

        ax.tick_params(axis='both', which='major', labelsize=11)

        # add subplot title (time)
        history_dt = GriddedOutput.historydelta_m/60 # hours
        time =  ((t-1)*history_dt)
        ax.set_title(f'$t = {time:3.0f}$ h', fontsize=11.5)

    if kwargs.get("savefig"):
        plt.savefig(f'speciated-mass-frac-four-panel-{scenario}-z{z_idx}.pdf', format='pdf', bbox_inches='tight')