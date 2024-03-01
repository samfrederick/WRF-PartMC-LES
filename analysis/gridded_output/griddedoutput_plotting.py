import matplotlib.pyplot as plt
from cycler import cycler
from load_griddedoutput_datasets import *
plt.rcParams['axes.prop_cycle'] = cycler('color', plt.get_cmap('tab20').colors)

def plotNumberDist(number_dist_arr, plot_type='hist', **kwargs):
    if not kwargs.get('ax'):
        fig, ax = plt.subplots(1, 1, figsize=(8,6))
    else:
        ax = kwargs.get('ax')

    if plot_type == 'line':
        c = kwargs.get('color', 'k')
        label = kwargs.get('label', None)
        lw = kwargs.get('lw', 1)
        ls = kwargs.get('ls', '-')
        ax.plot(bin_geocenter, number_dist_arr, c=c, label=label, lw=lw, ls=ls)
    elif plot_type == 'hist':
        bin_width = bin_edges[1:] - bin_edges[:-1] 
        ax.bar(x=bin_geocenter, height=number_dist_arr, width=bin_width, alpha=.6)
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

    if binned_species_mass_arr.shape[1] == n_bins:
        binned_species_mass_arr = binned_species_mass_arr[:, :-1]

    bin_total_mass = np.zeros(((n_bins-1),))
    fig, ax = plt.subplots(1, 1)
    plt.rcParams['axes.prop_cycle'] = cycler('color', plt.get_cmap('tab20').colors)

    # Plot the total speciated mass distribution
    for i, species in enumerate(aero_species):
        species_mass_dist = binned_species_mass_arr[i, :]#/bin_logwidth # moved division by logwidth to calculation of mass distrib
        bin_total_mass += species_mass_dist
        #print(bin_total_mass[50])
    #ax.plot(bin_edges[:-1], bin_total_mass, label=species)
        if i == 0:
            y_lower = ax.get_ylim()[0]
            y_lower = np.array((n_bins-1)*[y_lower])
        else:
            y_lower = bin_total_mass - species_mass_dist
        ax.fill_between(x=bin_geocenter, y1=y_lower, y2=bin_total_mass, label=species)
    
    plt.xscale('log')
    plt.yscale(kwargs.get('yscale', 'log'))
    plt.ylabel('Mass [kg$\cdot$m$^{-3}$]')
    plt.xlabel('Particle Diameter [m]')
    if kwargs.get('xlims'):
        plt.xlim(kwargs.get('xlims'))
    plt.legend(loc='center', bbox_to_anchor=(1.3, .5), ncol=2)
    if kwargs.get('time'):
        plt.text(1.3, .9, f't = {kwargs.get("time")} hrs', transform=ax.transAxes, horizontalalignment='center')
    if kwargs.get('savefig'):
        if kwargs.get('scenario'):
            plt.savefig(f'{kwargs.get("scenario")}_speciated_mass_dist_t{int(kwargs.get("time"))}.pdf', 
                        format='pdf', bbox_inches='tight')
        else:
            plt.savefig(f'speciated_mass_dist_t{int(kwargs.get("time"))}.pdf', format='pdf', bbox_inches='tight')

def plotSpeciatedMassFrac(binned_species_mass_arr, **kwargs):
    if binned_species_mass_arr.shape[1] == n_bins:
        binned_species_mass_arr = binned_species_mass_arr[:, :-1]

    total_mass_per_bin = (binned_species_mass_arr[:, :]#/bin_logwidth # moved division by logwidth to calculation of mass distrib
                          ).sum(axis=0)
    frac_total = np.zeros((n_bins-1,))
    fig, ax = plt.subplots(1, 1)
    # Plot the total speciated mass distribution
    for i, species in enumerate(aero_species):
        species_mass_dist = binned_species_mass_arr[i, :]#/bin_logwidth # moved division by logwidth to calculation of mass distrib
        species_frac = species_mass_dist/total_mass_per_bin
        frac_total += species_frac
        if i == 0:
            y_lower = 0
            y_lower = np.array((n_bins-1)*[y_lower])
        else:
            y_lower = frac_total - species_frac
        ax.fill_between(x=bin_geocenter, y1=y_lower, y2=frac_total, label=species)

    plt.xlim(bin_edges[0], bin_edges[-1])
    plt.ylim(0, 1)
    plt.xscale('log')
    plt.ylabel('Mass Fraction')
    plt.xlabel('Particle Diameter [m]')
    if kwargs.get('xlims'):
        plt.xlim(kwargs.get('xlims'))
    plt.legend(loc='center', bbox_to_anchor=(1.3, .5), ncol=2)
    if kwargs.get('time'):
        plt.text(1.3, .9, f't = {kwargs.get("time")} hrs', transform=ax.transAxes, horizontalalignment='center')
    if kwargs.get('savefig'):
        if kwargs.get('scenario'): 
            plt.savefig(f'{kwargs.get("scenario")}_speciated_mass_frac_t{int(kwargs.get("time"))}.pdf', 
                        format='pdf', bbox_inches='tight')
        else:
            plt.savefig(f'speciated_mass_frac_t{int(kwargs.get("time"))}.pdf', format='pdf', bbox_inches='tight')
