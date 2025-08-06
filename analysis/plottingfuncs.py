from helperfuncs import *
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mplcolors
import mpl_toolkits.axes_grid1.inset_locator as il
import matplotlib.patches as patches
import matplotlib.animation as animation
from matplotlib.ticker import FormatStrFormatter
from math import ceil
from scipy.ndimage.filters import gaussian_filter
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

def _getXTickTimes(xtick_units, xtick_delta, shift_tickloc=False):
    history_delta = Archive.historydelta_m
    total_duration = (Archive.n_times-1)*history_delta

    if xtick_units != 'm' and xtick_units != 'h':
        raise ValueError(f'Invalid label for xtick_units: {xtick_units}')
    elif xtick_units == 'h':
        total_duration = total_duration/60
        history_delta = history_delta/60

    n_xticks = int(1 + total_duration / xtick_delta)
    xtick_nudge = 0.5
    if not shift_tickloc:
        xtick_nudge = 0
    xticks = xtick_nudge + np.linspace(0, Archive.n_times-1, n_xticks) # shift by 1/2 to be cell-centered
    xtick_labels = history_delta*np.linspace(0, Archive.n_times-1, n_xticks)
    xtick_labels = xtick_labels.astype(int)

    return xticks, xtick_labels

def _createSuperSatLabel(variable):
    supersat = int(variable.split('_')[1])/10
    label = f'$S = {supersat}\%$'
    return label

def _plotSignificance(rel_diff, ax, thres_n_std_dev=5, skipzero=False):

    # get relative difference values during the first hour when the scenario matches up with the basecase
    # to get a baseline for how much variance should be expected due to statistical noise
    n_times_before_1hr = ceil(60/Archive.historydelta_m)
    if n_times_before_1hr < Archive.n_times:
        baseline = rel_diff[:n_times_before_1hr, :].flatten()
        if skipzero:
            baseline = rel_diff[1:n_times_before_1hr, :].flatten()

    # use a threshold of n*sigma (default to 5 sigma)
    thres = thres_n_std_dev*baseline.std()

    zm = np.ma.masked_less(abs(rel_diff), thres)

    ax.pcolor(zm.T, hatch='.', alpha=0.)

def _plotContours(ZT_data, ax, **kwargs):
    contour_smoothing = kwargs.get('smooth_contours', False)
    smoothing_sigma = kwargs.get('contour_smoothing_sigma', 0.7)
    contour_linewidth = kwargs.get('contour_linewidth', 1.5)
    contour_label_fontsize = kwargs.get('contour_label_fontsize', 9)
    contourmax = kwargs.get('contour_max', 15)
    contourmin = kwargs.get('contour_min', -15)
    ncontours = kwargs.get('n_contours', 9)
    contour_levels = kwargs.get('contour_levels', [x for x in np.linspace(contourmin,contourmax, ncontours) if x!= 0])
    contour_norm = kwargs.get('contour_norm', mplcolors.Normalize(contourmin, contourmax))
    fmt = kwargs.get('contour_fmt', matplotlib.ticker.ScalarFormatter())
    fmt.create_dummy_axis()
    # for whatever reason, the bounds for the countours are one cell off, so the contours on the 
    # rhs of plots was getting cut off. Shifting values to the right (one timestep) and copy paste
    # entries for the intiial condition into a new array with dimension (nlevels, ntimes+1)
    zt1 = ZT_data.T
    shift=kwargs.get('shift', 1)
    zt2 = np.zeros((zt1.shape[0], zt1.shape[1]+shift))
    zt2[:, shift:] = zt1[:, :]
    for i in range(shift):
        zt2[:, i] = zt1[:, i]

    if contour_smoothing:
        zt2 = gaussian_filter(zt2, smoothing_sigma)

    #CS = ax.contour(ZT_data.T, levels=contour_levels, cmap=kwargs.get('contour_cmap', 'gist_gray_r'))
    CS = ax.contour(zt2, levels=contour_levels, cmap=kwargs.get('contour_cmap', 'gist_gray_r'), 
                    norm=contour_norm, linewidths=contour_linewidth)
    ax.clabel(CS, inline=True, fontsize=contour_label_fontsize, fmt=fmt)

def plotNSH(scenario, variable, vmin=None, vmax=None, lognorm=False, **kwargs):
    
    if variable not in Archive.nsh_dict[scenario]:
        print(f'{variable} not in NSH dictionary for {scenario}, calculating')
        nsh_array = calculateNSHTimeSlice(scenario, variable)
    else:
        nsh_array = Archive.nsh_dict[scenario][variable]

    fig, ax  = plt.subplots(1,1, figsize=(12,5))
    if lognorm:
        norm = mplcolors.LogNorm(vmin, vmax)
    elif ((vmin != None) and (vmax != None)):
        norm = mplcolors.Normalize(vmin, vmax)
    else:
        norm = None
    cs = ax.pcolormesh(nsh_array.T, norm=norm,edgecolor='face')
    cbar = fig.colorbar(cs, label='NSH')

    ax.set_ylabel('z (km)', fontsize=12)

    # Set x-axis ticks and label
    xtick_units = kwargs.get('xtick_units', 'm') 
    xtick_delta = kwargs.get('xtick_delta_t', 30)
    xticks, xtick_labels = _getXTickTimes(xtick_units, xtick_delta, shift_tickloc=True)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xtick_labels)
    ax.set_xlabel(f'Time ({xtick_units})', fontsize=12)

    ax.set_yticks(np.arange(0, Archive.n_levels+1, 25))
    ax.set_yticklabels(np.linspace(0, 2, 5).round(2))
    ax.set_title(f'{scenario}, NSH ({variable})')

def plotMultipleNSH(scenario, variables, vmin=0, vmax=1, **kwargs):
    fig, axs  = plt.subplots(len(variables),1, figsize=kwargs.get('figsize', (6.5, 2.71)))

    axis_label_fontsize=10
    tick_label_fontsize=9

    if kwargs.get('lognorm', False):
        norm = mplcolors.LogNorm(vmin, vmax)
    elif ((vmin != None) and (vmax != None)):
        norm = mplcolors.Normalize(vmin, vmax)
    else:
        norm = None

    for i, (ax, variable) in enumerate(zip(axs.flatten(), variables)):
        variable_fmt = variable
        if variable in Archive.aero_vars:
            variable_fmt = Archive.aerosol_fmt_map[variable]
        elif variable in Archive.gas_vars:
            variable_fmt = Archive.gas_fmt_map[variable]

        if variable not in Archive.nsh_dict[scenario]:
            print(f'{variable} not in NSH dictionary for {scenario}, calculating')
            nsh_array = calculateNSHTimeSlice(scenario, variable)
        else:
            nsh_array = Archive.nsh_dict[scenario][variable]

        cs = ax.pcolormesh(nsh_array.T, norm=norm,edgecolor='face')


        ax.set_ylabel('z (km)', fontsize=axis_label_fontsize)

        # Set x-axis ticks and label
        xtick_units = kwargs.get('xtick_units', 'm') 
        xtick_delta = kwargs.get('xtick_delta_t', 30)
        xticks, xtick_labels = _getXTickTimes(xtick_units, xtick_delta, shift_tickloc=True)
        if i != len(variables)-1:
            ax.set_xticks([])
        else:
            ax.set_xticks(xticks)
            ax.set_xticklabels(xtick_labels, fontsize=tick_label_fontsize)
            ax.set_xlabel(f'Time ({xtick_units})', fontsize=axis_label_fontsize)

        ax.set_yticks(np.arange(0, Archive.n_levels+1, 25))
        ax.set_yticklabels(np.linspace(0, 2, 5).round(2), fontsize=tick_label_fontsize)
        #ax.set_title(f'$SH$, {variable_fmt}')
        ax.text(-.2, 0.5, f'{variable_fmt}', ha='center', transform=ax.transAxes, fontsize=axis_label_fontsize)

        if nsh_array.max() > 0.5:
            contourmax = 1
            contourmin = 0
            ncontours = 6
            kwargs['contour_max'] = contourmax
            kwargs['contour_min'] = contourmin
            kwargs['n_contours'] = ncontours
        elif (nsh_array.max() < 0.5) and (nsh_array.max() > 0.25):
            contourmax = .5
            contourmin = 0
            ncontours = 6
            kwargs['contour_max'] = contourmax
            kwargs['contour_min'] = contourmin
            kwargs['n_contours'] = ncontours
        elif (nsh_array.max() < 0.25) and (nsh_array.max() > 0.1):
            contourmax = .25
            contourmin = 0
            ncontours = 6
            kwargs['contour_max'] = contourmax
            kwargs['contour_min'] = contourmin
            kwargs['n_contours'] = ncontours
        else:
            contour_levels = np.logspace(-3, 0, 4)
            kwargs['contour_levels'] = contour_levels
        
        _plotContours(nsh_array, ax, **kwargs)
    
    fig.subplots_adjust(right=0.82)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(cs, cax=cbar_ax, label='$SH$')
    cbar.set_label('$SH$', fontsize=10)

    if kwargs.get("savefig"):
        plt.savefig(f'height-time-nsh-multivar-{scenario}.pdf', format='pdf', bbox_inches='tight')

def plotZT(scenario, variable, vmin=None, vmax=None, lognorm=False, **kwargs):
    
    fig, ax  = plt.subplots(1,1, figsize=(12,4.5))
    if lognorm:
        norm = mplcolors.LogNorm(vmin, vmax)
    elif (vmin and vmax):
        norm = mplcolors.Normalize(vmin, vmax)

    variable_fmt = variable
    if variable in Archive.aero_vars:
        if variable.startswith('ccn') or 'NUM_CONC' in variable:
            # number concentration to number per kg of air
            var_array = calculateVarZT(scenario, variable, convert_mixing_ratio=False)
            var_units = '# kg$^{-1}$'
        elif variable.startswith('pmc'):
            # convert to mixing ratio
            var_array = calculateVarZT(scenario, variable, convert_mixing_ratio=True)
            var_array = 1e9*var_array  # convert to ppbv
            var_units = 'Mixing Ratio (ppbv)'
        else:
            var_array = calculateVarZT(scenario, variable, convert_mixing_ratio=False)
            var_units = 'concentration'
        variable_fmt = Archive.aerosol_fmt_map[variable]
    elif variable in Archive.gas_vars:
        var_array = calculateVarZT(scenario, variable, convert_mixing_ratio=False)
        var_array = 1000*var_array # convert ppmv to ppbv
        var_units = 'Mixing Ratio (ppbv)'
        variable_fmt = Archive.gas_fmt_map[variable]
    elif variable in Archive.wrf_vars:
        var_array = calculateVarZT(scenario, variable, convert_mixing_ratio=False)
        var_units = 'concentration'

    cs = ax.pcolormesh(var_array.T, norm=norm, edgecolor='face')

    if kwargs.get("plot_contours"):
        _plotContours(var_array, ax, **kwargs)

    title = kwargs.get('title', f'{scenario}, {variable_fmt}')
    cbar_title = kwargs.get('cbar_title', f'{var_units}')
    cbar = fig.colorbar(cs)
    cbar.set_label(label=cbar_title, fontsize=13)
    cbar.ax.tick_params(labelsize=12)

    # Set x-axis ticks and label
    xtick_units = kwargs.get('xtick_units', 'm') 
    xtick_delta = kwargs.get('xtick_delta_t', 30)
    xticks, xtick_labels = _getXTickTimes(xtick_units, xtick_delta, shift_tickloc=True)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xtick_labels)
    ax.set_xlabel(f'Time ({xtick_units})', fontsize=14)

    # Set y-axis ticks and label
    ax.set_yticks(np.arange(0, Archive.n_levels+1, 25))
    ax.set_yticklabels(np.linspace(0, 2, 5).round(2))
    ax.set_ylabel('z (km)', fontsize=14)
    ax.tick_params(axis='both', which='major', labelsize=13)

    ax.set_title(title, fontsize=14)

    if kwargs.get("savefig"):
        plt.savefig(f'height-time-{variable}-{scenario}.pdf', format='pdf', bbox_inches='tight')

def plotFourPanelZT(scenarios, variable, vmin=None, vmax=None, lognorm=False, **kwargs):
    
    fig, axs  = plt.subplots(2,2, figsize=(6.5,4.5))
    plt.subplots_adjust(hspace=.3, wspace=.1)

    if lognorm:
        norm = mplcolors.LogNorm(vmin, vmax)
    elif (vmin and vmax):
        norm = mplcolors.Normalize(vmin, vmax)

    if len(scenarios) != 4:
        raise AttributeError("Number of scenarios must be four")
    
    labels= Archive.getScenarioGeneralLabels()
    
    for i, (ax, scenario) in enumerate(zip(axs.flatten(), scenarios)):

        variable_fmt = variable
        if variable in Archive.aero_vars:
            if variable.startswith('ccn') or 'NUM_CONC' in variable:
                # number concentration to number per kg of air
                var_array = calculateVarZT(scenario, variable, convert_mixing_ratio=False)
                var_units = '# kg$^{-1}$'
            elif variable.startswith('pmc'):
                # convert to mixing ratio
                var_array = calculateVarZT(scenario, variable, convert_mixing_ratio=True)
                var_array = 1e9*var_array  # convert to ppbv
                var_units = 'Mixing Ratio (ppbv)'
            else:
                var_array = calculateVarZT(scenario, variable, convert_mixing_ratio=False)
                var_units = 'concentration'
            variable_fmt = Archive.aerosol_fmt_map[variable]
        elif variable in Archive.gas_vars:
            var_array = calculateVarZT(scenario, variable, convert_mixing_ratio=False)
            var_array = 1000*var_array # convert ppmv to ppbv
            var_units = 'Mixing Ratio (ppbv)'
            variable_fmt = Archive.gas_fmt_map[variable]
        elif variable in Archive.wrf_vars:
            var_array = calculateVarZT(scenario, variable, convert_mixing_ratio=False)
            var_units = 'concentration'

        cs = ax.pcolormesh(var_array.T, norm=norm, edgecolor='face')

        if i==0:
            cbar = fig.colorbar(cs, ax=axs, orientation='horizontal', fraction=0.05, pad=0.15)
            cbar_title = kwargs.get('cbar_title', f'{var_units}')
            cbar.set_label(label=cbar_title, fontsize=11)
            cbar.ax.tick_params(labelsize=11)

        if kwargs.get("plot_contours"):
            _plotContours(var_array, ax, **kwargs)

        # Set x-axis ticks and labels   
        xtick_units = kwargs.get('xtick_units', 'm') 
        xtick_delta = kwargs.get('xtick_delta_t', 30)
        xticks, xtick_labels = _getXTickTimes(xtick_units, xtick_delta, shift_tickloc=True)
        ax.set_xticks(xticks)
        if i < 2:
            ax.set_xticklabels([])
        else:
            ax.set_xticklabels(xtick_labels)
            ax.set_xlabel(f'Time ({xtick_units})', fontsize=11)

        # Set y-axis ticks and labels
        ax.set_yticks(np.arange(0, Archive.n_levels+1, 25))
        if i in [1,3]:
            ax.set_yticklabels([])
        else:
            # Set y-axis ticks and label
            ax.set_yticklabels(np.linspace(0, 2, 5).round(2))
            ax.set_ylabel('z (km)', fontsize=11)

        ax.tick_params(axis='both', which='major', labelsize=11)

        title = labels[scenario]
        ax.set_title(title, fontsize=11.5)

    if kwargs.get("savefig"):
        plt.savefig(f'height-time-{variable}-four-scenarios.pdf', format='pdf', bbox_inches='tight')

def plotConcT(scenario, variable, zlevel, vmin=None, vmax=None, lognorm=False, convert_mixing_ratio=True, **kwargs):
    
    var_array = calculateVarZT(scenario, variable, convert_mixing_ratio)
    
    fig, ax  = plt.subplots(1,1, figsize=(12,5))
    if lognorm:
        norm = mplcolors.LogNorm(vmin, vmax)
    else:
        norm = None
    cs = ax.pcolormesh(var_array.T, norm=norm,edgecolor='face')

    if convert_mixing_ratio:
        var_units = 'mixing ratio'
    else:
        var_units = 'concentration'
    
    title = kwargs.get('title', f'{scenario}, {variable}')
    cbar_title = kwargs.get('cbar_title', f'{variable}, {var_units}')
    cbar = fig.colorbar(cs, label=cbar_title)

    # Set x-axis ticks and label
    xtick_units = kwargs.get('xtick_units', 'm') 
    xtick_delta = kwargs.get('xtick_delta_t', 30)
    xticks, xtick_labels = _getXTickTimes(xtick_units, xtick_delta)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xtick_labels)
    ax.set_xlabel(f'Time ({xtick_units})', fontsize=12)

    ax.set_ylabel('z (km)', fontsize=12)
    ax.set_yticks(np.arange(0, Archive.n_levels+1, 25))
    ax.set_yticklabels(np.linspace(0, 2, 5).round(2))
    ax.set_title(title)

def plotScenariosVarsLevelConc(scenarios, variables, zlevel, convert_mixing_ratio, **kwargs):
    #var_array = calculateVarZT(scenario, variable, convert_mixing_ratio)
    n_vars = len(variables)
    if not kwargs.get('ax'):
        fig, axs  = plt.subplots(n_vars,1, figsize=(kwargs.get('subplot_width', 12),n_vars*kwargs.get('subplot_height', 3)))
    else:
        axs = kwargs.get('ax')
    if kwargs.get('general_scenario_label'):
        labels= Archive.getScenarioGeneralLabels()
        colors = Archive.getScenarioGeneralColors()   

    if not isinstance(axs, np.ndarray):
        axs = np.array([axs])
    for ax, variable in zip(axs, variables):
        for scenario in scenarios:
            times = np.arange(Archive.n_times)
            var_array = np.zeros((Archive.n_times))
            for itime in times:
                if convert_mixing_ratio:
                    inverse_airdens = Archive.aero_data[scenario]['ALT'][itime, zlevel, :, :]
                    level_array = inverse_airdens*Archive.aero_data[scenario][variable][itime, zlevel, :, :]
                else:
                    level_array = Archive.aero_data[scenario][variable][itime, zlevel, :, :]
                var_array[itime] = level_array.mean()

            if kwargs.get('general_scenario_label'):
                label = labels[scenario]
                color = colors[scenario]
                ax.plot(times, var_array, label = label, c=color)
            else:
                label = scenario
                ax.plot(times, var_array, label = label)
            

        if convert_mixing_ratio:
            var_units = 'Mixing Ratio'
        else:
            var_units = 'Concentration'

        #cbar = fig.colorbar(cs, label=f'{variable} {var_units}')
        if kwargs.get('plot_legend', True):
            ax.legend()
        # Set x-axis ticks and label
        xtick_units = kwargs.get('xtick_units', 'm') 
        xtick_delta = kwargs.get('xtick_delta_t', 30)
        xticks, xtick_labels = _getXTickTimes(xtick_units, xtick_delta)
        ax.set_xticks(xticks)
        ax.set_xticklabels(xtick_labels)
        ax.set_xlabel(f'Time ({xtick_units})', fontsize=12)

        ax.set_ylabel(f'{var_units}', fontsize=12)
        #ax.set_yticks(np.arange(0, Archive.n_levels+1, 25))
        #ax.set_yticklabels(np.linspace(0, 2, 5).round(2))
        ax.set_title(f'{variable}')

        if kwargs.get('grid', True):
            ax.set_xlim(0, (Archive.n_times-1))
            ax.grid(which = "major", linewidth = 1, axis='y', ls="dotted", dashes=(.5,6), c='#414141', alpha=.5)
            ax.grid(which = "minor", linewidth = 1, axis='y', ls="dotted", dashes=(.5,6), c='white')
            ax.grid(which = "minor", linewidth = 1, axis='x', ls="dotted", dashes=(.5,6), c='#414141')
            ax.grid(which = "major", linewidth = 1, axis='x', ls="dotted", dashes=(.5,6), c='#414141')
            ax.tick_params(axis='both', labelsize=13, which='major', direction='in', top=True, right=True, bottom=True, left=True)
            ax.tick_params(axis='both', which='minor',direction='in',top=True, right=True, bottom=True, left=True)

    plt.suptitle(f'Z-level: {zlevel}')
    plt.tight_layout()

def plotScenariosVarsVerticalProfile(scenarios, variables, time, **kwargs):
    #var_array = calculateVarZT(scenario, variable, convert_mixing_ratio)
    global_fontsize = kwargs.get('global_fontsize', 13)
    if (len(variables) < 4):
        fig_xsize = kwargs.get('fig_xsize', len(variables)*4)
        fig_ysize = kwargs.get('fig_ysize', 4.5)
        fig, axs  = plt.subplots(1,len(variables), figsize=(fig_xsize,fig_ysize))

    if (len(variables) == 4):
        fig_xsize = kwargs.get('fig_xsize', 2*4)
        fig_ysize = kwargs.get('fig_ysize', 8.5)
        fig, axs  = plt.subplots(2, 2, figsize=(fig_xsize,fig_ysize))

    if (len(variables) > 4) and (len(variables) < 7):
        fig_xsize = kwargs.get('fig_xsize', 3*4)
        fig_ysize = kwargs.get('fig_ysize', 8.5)
        fig, axs  = plt.subplots(2, 3, figsize=(fig_xsize,fig_ysize))
        if len(variables)<6:
            axes_to_remove = 6-len(variables)
            for i in range(axes_to_remove):
                i+=1
                axs[-1*i].remove()

    if kwargs.get('general_scenario_label'):
        labels= Archive.getScenarioGeneralLabels()
    if kwargs.get('use_standard_scenario_colors'):
        colors = Archive.getScenarioColors()
    else:
        set1 = plt.get_cmap(kwargs.get('cmap', 'Set1'))
        num_colors = set1.N # number of colors in the set1 colormap
        colors = {}
        # Extract colors and convert to hex
        for i, scenario in enumerate(scenarios): # alternatively loop over range num_colors
            rgba = set1(i)  # Get RGBA values
            hex_color = mplcolors.to_hex(rgba)  # Convert RGBA to hex
            colors[scenario] = hex_color

    vars_to_quantile = kwargs.get('vars_to_quantile', {})

    for i, (ax, variable) in enumerate(zip(axs.flatten(), variables)):
        write_quantile=True
        for scenario in scenarios:
            #times = np.arange(Archive.n_times)
            var_array = np.zeros((Archive.n_times))
            variable_fmt = variable
            if variable in Archive.aero_vars:
                if variable.startswith('ccn') or 'NUM_CONC' in variable:
                    # number concentration to number per kg of air
                    inverse_airdens = Archive.aero_data[scenario]['ALT'][time, :, :, :]
                    array = inverse_airdens*Archive.aero_data[scenario][variable][time, :, :, :]
                    var_units = '# kg$^{-1}$'
                    if kwargs.get('unit_prefactor'):
                        unit_prefactor = kwargs.get('unit_prefactor')
                        array = unit_prefactor*inverse_airdens*Archive.aero_data[scenario][variable][time, :, :, :]
                        unit_prefactor_str = f'{1/unit_prefactor:1.0e}'
                        components = unit_prefactor_str.split('e')
                        signif = components[0]
                        expon = components[1].replace('+', '').lstrip("0")
                        fmt_unit_prefac = '$' + signif + '\\times' + '10^{' + expon + '}' + '$'

                        var_units = '(# kg$^{-1}$)$/$' + fmt_unit_prefac

                elif variable.startswith('pmc'):
                    # convert to mixing ratio
                    inverse_airdens = Archive.aero_data[scenario]['ALT'][time, :, :, :]
                    array = 1e9*inverse_airdens*Archive.aero_data[scenario][variable][time, :, :, :]
                    var_units = 'Mixing Ratio (ppbv)'

                else:
                    array = Archive.aero_data[scenario][variable][time, :, :, :]
                    var_units = ''
                variable_fmt = Archive.aerosol_fmt_map[variable]
            elif variable in Archive.gas_vars:
                if variable == 'oh':
                    scaling_factor = 1e6 # scale to pptv
                    var_units = 'Mixing Ratio (pptv)'
                else:
                    scaling_factor = 1e3
                    var_units = 'Mixing Ratio (ppbv)'
                array = scaling_factor*Archive.aero_data[scenario][variable][time, :, :, :] # convert ppmv to ppbv
                
                variable_fmt = Archive.gas_fmt_map[variable]
            elif variable in Archive.wrf_vars:
                array = Archive.aero_data[scenario][variable][time, :, :, :]
                var_units = ''
            else:
                raise AttributeError(f'Variable {variable} not recognized')
            
            if variable in vars_to_quantile.keys():
                quantile = vars_to_quantile[variable]
                array_mean = np.quantile(array, quantile, axis=(1,2))
                if write_quantile:
                    ax.text(.5, .87, f'{quantile} quantile', transform=ax.transAxes, fontsize=global_fontsize, horizontalalignment='center', 
                            bbox=dict(facecolor='white', edgecolor='black', alpha=.8))
                    write_quantile = False
            else:
                array_mean = array.mean(axis=(1,2))

            if kwargs.get('general_scenario_label'):
                label = labels[scenario]
            else:
                label = scenario
            
            color = colors[scenario]
            ls = '-'
            if 'no-nh4' in scenario:
                ls = '--'
            ax.plot(array_mean, np.arange(100), label=label, c=color, ls=ls)
            #if kwargs.get('plot_std'):
            #    array_std = array.std(axis=(1,2))
            #    ax.fill_betweenx(np.arange(100), array_mean-array_std, array_mean+array_std, alpha=.4)

        #cbar = fig.colorbar(cs, label=f'{variable} {var_units}')
        legend_cols = kwargs.get('legend_ncols', len(scenarios))
        if len(variables) > 1:
            if i == 0:
                fig.legend(fontsize=global_fontsize, ncol=legend_cols, loc='center', bbox_to_anchor=(.5,-.05))
        else:
            ax.legend(fontsize=global_fontsize)
        ax.set_xlabel(f'{var_units}', fontsize=global_fontsize)
        ax.set_ylabel('z (km)', fontsize=global_fontsize)
        ax.set_yticks(np.arange(0, Archive.n_levels+1, 25))
        ax.set_yticklabels(np.linspace(0, 2, 5).round(2))
        ax.set_title(f'{variable_fmt}', fontsize=global_fontsize)

        ax.set_ylim(0, Archive.n_levels)
        if kwargs.get('grid', True):
            lw = kwargs.get('grid_linewidth', 1)
            ax.grid(which = "major", linewidth = lw, axis='y', ls="dotted", dashes=(.5,6), c='#414141', alpha=.5)
            ax.grid(which = "minor", linewidth = lw, axis='y', ls="dotted", dashes=(.5,6), c='white')
            ax.grid(which = "minor", linewidth = lw, axis='x', ls="dotted", dashes=(.5,6), c='#414141')
            ax.grid(which = "major", linewidth = lw, axis='x', ls="dotted", dashes=(.5,6), c='#414141')
            ax.tick_params(axis='both', labelsize=global_fontsize, which='major', direction='in', top=True, right=True, bottom=True, left=True)
            ax.tick_params(axis='both', which='minor',direction='in',top=True, right=True, bottom=True, left=True)

    delta_t = Archive.historydelta_m
    plt.suptitle(kwargs.get('title', f'Time: {delta_t*time} m'))
    plt.tight_layout()

    if kwargs.get("savefig"):
        filename = kwargs.get('filename', f'vertical-profiles-time{time}.pdf')
        plt.savefig(filename, format='pdf', bbox_inches='tight')

def plotNSHPercentDiff(scenario, variable, vmin=None, vmax=None, **kwargs):
    
    if variable not in Archive.nsh_dict[scenario]:
        print(f'{variable} not in NSH dictionary for {scenario}, calculating')
        nsh_array_scenario = calculateNSHTimeSlice(scenario, variable)
    else:
        nsh_array_scenario = Archive.nsh_dict[scenario][variable]

    if variable not in Archive.nsh_dict['uniform-basecase']:
        print(f'{variable} not in NSH dictionary for uniform-basecase, calculating')
        nsh_array_basecase = calculateNSHTimeSlice('uniform-basecase', variable)
    else:
        nsh_array_basecase = Archive.nsh_dict['uniform-basecase'][variable]
    
    rel_diff = 100*(nsh_array_scenario - nsh_array_basecase)/nsh_array_basecase
    fig, ax  = plt.subplots(1,1, figsize=(12,5))
    cs = ax.pcolormesh(rel_diff.T, cmap=plt.cm.coolwarm, vmin=vmin, vmax=vmax,edgecolor='face')
    cbar = fig.colorbar(cs, label='NSH percent difference')

    # Set x-axis ticks and label
    xtick_units = kwargs.get('xtick_units', 'm') 
    xtick_delta = kwargs.get('xtick_delta_t', 30)
    xticks, xtick_labels = _getXTickTimes(xtick_units, xtick_delta, shift_tickloc=True)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xtick_labels)
    ax.set_xlabel(f'Time ({xtick_units})', fontsize=12)

    ax.set_ylabel('z (km)', fontsize=12)
    ax.set_yticks(np.arange(0, Archive.n_levels+1, 25))
    ax.set_yticklabels(np.linspace(0, 2, 5).round(2))
    ax.set_title(f'{scenario} NSH ({variable}) percent difference')

def plotVarPercentDiff(scenario, variable, vmin=None, vmax=None, convert_mixing_ratio=False, skip_t0=False, **kwargs):
    
    rel_diff = calculateVarPercentDiff(scenario, variable, convert_mixing_ratio, skip_t0)
    fig, ax  = plt.subplots(1,1, figsize=(12,4.5))
    cs = ax.pcolormesh(rel_diff.T, cmap=plt.cm.coolwarm, vmin=vmin, vmax=vmax,edgecolor='face')
    cbar = fig.colorbar(cs)
    cbar_title = kwargs.get('cbar_title', '% difference')   
    cbar.set_label(label=cbar_title, fontsize=13)
    cbar.ax.tick_params(labelsize=12)

    if kwargs.get("plot_contours"):
        _plotContours(rel_diff, ax, **kwargs)

    variable_fmt = variable
    if variable in Archive.aero_vars:
        variable_fmt = Archive.aerosol_fmt_map[variable]
    if variable in Archive.gas_vars:
        variable_fmt = Archive.gas_fmt_map[variable]

    # plot significance hatching
    if kwargs.get('plot_significance'):
        _plotSignificance(rel_diff, ax, thres_n_std_dev=kwargs.get('signif_thres_n_std_dev', 5), skipzero=skip_t0)

    # Set x-axis ticks and label
    xtick_units = kwargs.get('xtick_units', 'm') 
    xtick_delta = kwargs.get('xtick_delta_t', 30)
    xticks, xtick_labels = _getXTickTimes(xtick_units, xtick_delta, shift_tickloc=True)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xtick_labels)
    ax.set_xlabel(f'Time ({xtick_units})', fontsize=14)

    # Set y-axis ticks and label
    ax.set_yticks(np.arange(0, Archive.n_levels+1, 25))
    ax.set_yticklabels(np.linspace(0, 2, 5).round(2))
    ax.set_ylabel('z (km)', fontsize=14)
    ax.tick_params(axis='both', which='major', labelsize=13)

    title = kwargs.get('title', f'{scenario} {variable_fmt} % difference ')
    ax.set_title(title, fontsize=14)

    if kwargs.get("savefig"):
        plt.savefig(f'height-time-pdiff-{variable}-{scenario}.pdf', format='pdf', bbox_inches='tight')

def plotFourPanelVarPercentDiff(scenarios, variable, vmin=None, vmax=None, convert_mixing_ratio=False, skip_t0=False, **kwargs):
    
    fig, axs  = plt.subplots(2,2, figsize=(6.5,4.5))
    plt.subplots_adjust(hspace=.3, wspace=.1)

    if len(scenarios) != 4:
        raise AttributeError("Number of scenarios must be four")
    
    labels= Archive.getScenarioGeneralLabels()

    for i, (ax, scenario) in enumerate(zip(axs.flatten(), scenarios)):
        rel_diff = calculateVarPercentDiff(scenario, variable, convert_mixing_ratio, skip_t0)
        
        cs = ax.pcolormesh(rel_diff.T, cmap=plt.cm.coolwarm, vmin=vmin, vmax=vmax,edgecolor='face')

        if i==0:
            cbar = fig.colorbar(cs, ax=axs, orientation='horizontal', fraction=0.05, pad=0.15
                                )

            cbar_title = kwargs.get('cbar_title', '% difference')   
            cbar.set_label(label=cbar_title, fontsize=11)
            cbar.ax.tick_params(labelsize=11)

        if kwargs.get("plot_contours"):
            _plotContours(rel_diff, ax, **kwargs)

        variable_fmt = variable
        if variable in Archive.aero_vars:
            variable_fmt = Archive.aerosol_fmt_map[variable]
        if variable in Archive.gas_vars:
            variable_fmt = Archive.gas_fmt_map[variable]

        # plot significance hatching
        if kwargs.get('plot_significance'):
            _plotSignificance(rel_diff, ax, thres_n_std_dev=kwargs.get('signif_thres_n_std_dev', 5), skipzero=skip_t0)

        # Set x-axis ticks and labels   
        xtick_units = kwargs.get('xtick_units', 'm') 
        xtick_delta = kwargs.get('xtick_delta_t', 30)
        xticks, xtick_labels = _getXTickTimes(xtick_units, xtick_delta, shift_tickloc=True)
        ax.set_xticks(xticks)
        if i < 2:
            ax.set_xticklabels([])
        else:
            ax.set_xticklabels(xtick_labels)
            ax.set_xlabel(f'Time ({xtick_units})', fontsize=11)

        # Set y-axis ticks and labels
        ax.set_yticks(np.arange(0, Archive.n_levels+1, 25))
        if i in [1,3]:
            ax.set_yticklabels([])
        else:
            # Set y-axis ticks and label
            ax.set_yticklabels(np.linspace(0, 2, 5).round(2))
            ax.set_ylabel('z (km)', fontsize=11)

        ax.tick_params(axis='both', which='major', labelsize=11)

        title = labels[scenario]
        ax.set_title(title, fontsize=11.5)

    if kwargs.get("savefig"):
        plt.savefig(f'height-time-{variable}-pdiff-four-scenarios.pdf', format='pdf', bbox_inches='tight')

def plotCCNPercentDiff(scenario, vmin=None, vmax=None, convert_mixing_ratio=False, skip_t0=False, **kwargs):
    
    fig, axs  = plt.subplots(2,2, figsize=(6.5,4.5))
    plt.subplots_adjust(hspace=.3, wspace=.1)
    ccn_vars = ['ccn_001', 'ccn_003', 'ccn_006', 'ccn_010']

    for i, (ax, variable) in enumerate(zip(axs.flatten(), ccn_vars)):
        rel_diff = calculateVarPercentDiff(scenario, variable, convert_mixing_ratio, skip_t0)
        
        cs = ax.pcolormesh(rel_diff.T, cmap=plt.cm.coolwarm, vmin=vmin, vmax=vmax,edgecolor='face')

        if i==0:
            cbar = fig.colorbar(cs, ax=axs, orientation='horizontal', fraction=0.05, pad=0.15
                                )

            cbar_title = kwargs.get('cbar_title', '% difference')   
            cbar.set_label(label=cbar_title, fontsize=11)
            cbar.ax.tick_params(labelsize=11)

        if kwargs.get("plot_contours"):
            _plotContours(rel_diff, ax, **kwargs)

        variable_fmt = variable
        if variable in Archive.aero_vars:
            variable_fmt = Archive.aerosol_fmt_map[variable]
        if variable in Archive.gas_vars:
            variable_fmt = Archive.gas_fmt_map[variable]

        # plot significance hatching
        if kwargs.get('plot_significance'):
            _plotSignificance(rel_diff, ax, thres_n_std_dev=kwargs.get('signif_thres_n_std_dev', 5), skipzero=skip_t0)

        # Set x-axis ticks and labels   
        xtick_units = kwargs.get('xtick_units', 'm') 
        xtick_delta = kwargs.get('xtick_delta_t', 30)
        xticks, xtick_labels = _getXTickTimes(xtick_units, xtick_delta, shift_tickloc=True)
        ax.set_xticks(xticks)
        if i < 2:
            ax.set_xticklabels([])
        else:
            ax.set_xticklabels(xtick_labels)
            ax.set_xlabel(f'Time ({xtick_units})', fontsize=11)

        # Set y-axis ticks and labels
        ax.set_yticks(np.arange(0, Archive.n_levels+1, 25))
        if i in [1,3]:
            ax.set_yticklabels([])
        else:
            # Set y-axis ticks and label
            ax.set_yticklabels(np.linspace(0, 2, 5).round(2))
            ax.set_ylabel('z (km)', fontsize=11)

        ax.tick_params(axis='both', which='major', labelsize=11)

        #title = kwargs.get('title', f'{scenario} {variable_fmt} % difference ')
        ax.set_title(variable_fmt, fontsize=11.3)

    if kwargs.get("savefig"):
        plt.savefig(f'height-time-ccn-pdiff-{scenario}.pdf', format='pdf', bbox_inches='tight')


def plotMultiScenarioCCNPercentDiff(vmin=None, vmax=None, convert_mixing_ratio=False, skip_t0=False, **kwargs):

    fig_xsize = kwargs.get('fig_xsize', 7)
    fig_ysize = kwargs.get('fig_ysize', 3)
    fig, axs  = plt.subplots(4,3, figsize=(fig_xsize, fig_ysize), sharex=True, sharey=True, layout='constrained')

    global_fontsize = kwargs.get('global_fontsize', 9)

    general_scenario_labels = Archive.getScenarioGeneralLabels()

    vars = ['ccn_001', 'ccn_003', 'ccn_006', 'ccn_010']
    scenarios = ['fx1fy0', 'road-10x', 'point-source-1x1']*4
    for i, (scenario, ax) in enumerate(zip(scenarios, axs.flatten())):
        if i % 3 == 0:
            variable = vars[i//3]
        
        rel_diff = calculateVarPercentDiff(scenario, variable, convert_mixing_ratio, skip_t0)

        cs = ax.pcolormesh(rel_diff.T, cmap=plt.cm.coolwarm, vmin=vmin, vmax=vmax,edgecolor='face')

        if i==0:
            cbar = fig.colorbar(cs, ax=axs, orientation='horizontal', fraction=0.04, pad=0.05
                                )

            cbar_title = kwargs.get('cbar_title', '% difference')   
            cbar.set_label(label=cbar_title, fontsize=global_fontsize)
            cbar.ax.tick_params(labelsize=global_fontsize)

        if kwargs.get("plot_contours"):
            _plotContours(rel_diff, ax, **kwargs)

        variable_fmt = variable
        if variable in Archive.aero_vars:
            variable_fmt = Archive.aerosol_fmt_map[variable]
        if variable in Archive.gas_vars:
            variable_fmt = Archive.gas_fmt_map[variable]

        # plot significance hatching
        if kwargs.get('plot_significance'):
            _plotSignificance(rel_diff, ax, thres_n_std_dev=kwargs.get('signif_thres_n_std_dev', 5), skipzero=skip_t0)

        # Set x-axis ticks and labels   
        xtick_units = kwargs.get('xtick_units', 'm') 
        xtick_delta = kwargs.get('xtick_delta_t', 30)
        xticks, xtick_labels = _getXTickTimes(xtick_units, xtick_delta, shift_tickloc=True)
        ax.set_xticks(xticks)
        if i < 9:
            ax.set_xticklabels([])
        else:
            ax.set_xticklabels(xtick_labels)
            ax.set_xlabel(f'Time ({xtick_units})', fontsize=global_fontsize)

        # Set y-axis ticks and labels
        ax.set_yticks(np.arange(0, Archive.n_levels+1, 25))
        ax.set_yticklabels(np.linspace(0, 2, 5).round(2))
        if i%3 == 0:
            ax.set_ylabel('z (km)', fontsize=global_fontsize)

            # add text 
            label = _createSuperSatLabel(variable)
            ax.text(x=-.6, y=.5, s=label, transform=ax.transAxes)

        ax.tick_params(axis='both', which='major', labelsize=global_fontsize)

        if i < 3:
            ax.set_title(general_scenario_labels[scenario], fontsize=global_fontsize)

    if kwargs.get("savefig"):
        plt.savefig(f'height-time-ccn-pdiff-multi-scenario.pdf', format='pdf', bbox_inches='tight')


def plotVarBias(scenario, variable, vmin=None, vmax=None, convert_mixing_ratio=False, skip_t0=False, **kwargs):
    
    bias = calculateVarBias(scenario, variable, convert_mixing_ratio, skip_t0)
    fig, ax  = plt.subplots(1,1, figsize=(12,5))
    cs = ax.pcolormesh(bias.T, cmap=plt.cm.coolwarm, vmin=vmin, vmax=vmax,edgecolor='face')
    cbar = fig.colorbar(cs, label=f'{variable} bias')

    # Set x-axis ticks and label
    xtick_units = kwargs.get('xtick_units', 'm') 
    xtick_delta = kwargs.get('xtick_delta_t', 30)
    xticks, xtick_labels = _getXTickTimes(xtick_units, xtick_delta, shift_tickloc=True)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xtick_labels)
    ax.set_xlabel(f'Time ({xtick_units})', fontsize=12)

    ax.set_ylabel('z (km)', fontsize=12)
    ax.set_yticks(np.arange(0, Archive.n_levels+1, 25))
    ax.set_yticklabels(np.linspace(0, 2, 5).round(2))

    if convert_mixing_ratio:
        mixingratio_str = 'mixing ratio '
    else:
        mixingratio_str = ''
    ax.set_title(f'{scenario} {variable} {mixingratio_str}bias')

def plotNumberDist(scenario, i, j, k, **kwargs):
    scenario_aerodata = Archive.aero_data[scenario]
    scenario_distdata = Archive.aerodist_data[scenario]
    fig, ax = plt.subplots(1,2, figsize=(8,4))

    # Configurable keyword arguments
    times = kwargs.get('times', np.arange(0, Archive.n_times+1, int(Archive.n_times/(60/Archive.historydelta_m))))
    xlims = kwargs.get('xlims', (5e-9, 5e-6))
    ylims = kwargs.get('ylims', (0, 10000))
    numconctimeidx=kwargs.get('numconctimeidx', Archive.n_times-1)
    local_binning = kwargs.get('local_binning', None)
    if not isinstance(times, np.ndarray):
        times = np.array(times)

    cmap_name = kwargs.get('cmap')
    if not cmap_name:
        cmap = plt.cm.viridis
    else:
        cmap =plt.get_cmap(cmap_name)

    colors = cmap(np.linspace(.2, .9, times.size))
    for c, time in zip(colors, times):
        x_vals = []
        bin_vals = []
        bin_edges = scenario_distdata['BIN_EDGES'][:].data[0]#scenario_aerodata['BIN_EDGES'][:].data[0]
        bin_centers = scenario_distdata['BIN_CENTERS'][:].data[0]#scenario_aerodata['BIN_CENTERS'][:].data[0]
        bin_width = bin_edges[1:] - bin_edges[:-1]
        for bin_idx in range(100):
            bin_idx += 1 # 1 indexing 
            bin_data = scenario_distdata[f'num_a{str(bin_idx).zfill(3)}'][time, k, j, i].data.item()#/1e6
            if local_binning:
                bin_data = scenario_distdata[f'num_a{str(bin_idx).zfill(3)}'][time, k, j-local_binning:j+local_binning, i-local_binning:i+local_binning].data#/1e6
                bin_data = bin_data.mean()
            bin_vals.append(bin_data)
            x_vals.append(bin_idx)

        ax[0].plot(bin_centers, bin_vals, label=f't={Archive.historydelta_m*time}', c=c, lw=1.5)
        ax[0].set_xscale('log')

    ax[0].set_xlim(xlims[0], xlims[1])
    ax[0].set_ylim(ylims[0], ylims[1])
    ax[0].legend()
    ax[0].set_ylabel('Number concentration [m$^{-3}$]')
    ax[0].set_xlabel('Diameter [m]')
    ax[0].set_title(f'Number distribution')

    ax[1].pcolormesh(scenario_aerodata['TOT_NUM_CONC'][numconctimeidx, k, :, :])
    ax[1].plot([i], [j], marker='*', c='white', markeredgecolor='k', markeredgewidth=.5, markersize=15)
    ax[1].set_ylabel('j')
    ax[1].set_xlabel('i')
    ax[1].set_title(f'Total number conc. at t = {Archive.historydelta_m*numconctimeidx} m')

    if local_binning:
        ax[1].add_patch(
            patches.Rectangle(
                (i-local_binning, j-local_binning),
                2*local_binning,
                2*local_binning,
                fill=True,      # remove background
                color='white',
                alpha=.5
            ) ) 

    plt.suptitle(f'i={i}, j={j}, k={k}')
    plt.tight_layout()

    ax[0].set_box_aspect(1)
    ax[1].set_box_aspect(1)
    
def plotNumberDists(scenarios, i, j, k, time, **kwargs):
    
    fig, ax = plt.subplots(1,2, figsize=(8,4))

    # Configurable keyword arguments
    
    xlims = kwargs.get('xlims', (5e-9, 5e-6))
    ylims = kwargs.get('ylims', (0, 10000))
    numconctimeidx=kwargs.get('numconctimeidx', Archive.n_times-1)
    local_binning = kwargs.get('local_binning', None)
    #colors = plt.cm.viridis(np.linspace(.2, .9, times.size))

    for scenario in scenarios:
        scenario_aerodata = Archive.aero_data[scenario]
        scenario_distdata = Archive.aerodist_data[scenario]


        x_vals = []
        bin_vals = []
        bin_edges = scenario_distdata['BIN_EDGES'][:].data[0]#scenario_aerodata['BIN_EDGES'][:].data[0]
        bin_centers = scenario_distdata['BIN_CENTERS'][:].data[0]#scenario_aerodata['BIN_CENTERS'][:].data[0]
        bin_width = bin_edges[1:] - bin_edges[:-1]
        for bin_idx in range(100):
            bin_idx += 1 # 1 indexing 
            bin_data = scenario_distdata[f'num_a{str(bin_idx).zfill(3)}'][time, k, j, i].data.item()#/1e6
            if local_binning:
                bin_data = scenario_distdata[f'num_a{str(bin_idx).zfill(3)}'][time, k, j-local_binning:j+local_binning, i-local_binning:i+local_binning].data#/1e6
                bin_data = bin_data.mean()
            bin_vals.append(bin_data)
            x_vals.append(bin_idx)

        ax[0].plot(bin_centers, bin_vals, label=f'{scenario}', lw=1.5)
        ax[0].set_xscale('log')

    ax[0].set_xlim(xlims[0], xlims[1])
    ax[0].set_ylim(ylims[0], ylims[1])
    ax[0].legend()
    ax[0].set_ylabel('Number concentration [m$^{-3}$]')
    ax[0].set_xlabel('Diameter [m]')
    ax[0].set_title('Number distribution')

    """
    ax[1].pcolormesh(scenario_aerodata['TOT_NUM_CONC'][numconctimeidx, k, :, :])
    ax[1].plot([i], [j], marker='*', c='white', markeredgecolor='k', markeredgewidth=.5, markersize=15)
    ax[1].set_ylabel('j')
    ax[1].set_xlabel('i')
    ax[1].set_title(f'Total number conc. at t = {Archive.historydelta_m*numconctimeidx} m')

    if local_binning:
        ax[1].add_patch(
            patches.Rectangle(
                (i-local_binning, j-local_binning),
                2*local_binning,
                2*local_binning,
                fill=True,      # remove background
                color='white',
                alpha=.5
            ) ) 
    """
    plt.suptitle(f'i={i}, j={j}, k={k}, t={time}')

def plotDistand3DCrossSec(scenario, i, j, k, dist_type, **kwargs):
    # https://matplotlib.org/stable/gallery/mplot3d/mixed_subplots.html
    # https://matplotlib.org/stable/gallery/mplot3d/2dcollections3d.html

    scenario_aerodata = Archive.aero_data[scenario]
    scenario_distdata = Archive.aerodist_data[scenario]

    # Configurable keyword arguments
    times = kwargs.get('times', np.arange(0, Archive.n_times+1, int(Archive.n_times/(60/Archive.historydelta_m))))
    xlims = kwargs.get('xlims', (5e-9, 5e-6))
    if dist_type == 'num':
        ylims = kwargs.get('ylims', (1e7, 1e10))
        scaling_factor = 1 # keep as # m^-3
    if dist_type == 'mass':
        ylims = kwargs.get('ylims', (0, 1e-13))
        scaling_factor=1e9 # convert to micrograms per cubic meter
    yscale = kwargs.get('yscale', 'linear')
    savefig = kwargs.get('savefig', True)
    totconctimeidx=kwargs.get('totconctimeidx', Archive.n_times-1)
    lognorm = kwargs.get('lognorm', False)
    if not isinstance(times, np.ndarray):
        times = np.array(times)
    dist_cmap_name = kwargs.get('dist_cmap', 'viridis')
    dist_cmap = plt.get_cmap(dist_cmap_name)
    dist_cmap_normrange = kwargs.get('dist_cmap_norm_range', (0.2, 0.9))
    title_size = kwargs.get('title_size', 12)
    label_size = kwargs.get('label_size', 11)
    local_binning = kwargs.get('local_binning', None)
    legend_loc = kwargs.get('legend_loc', 'upper left')
    legend_ncol = kwargs.get('legend_ncol', len(times))

    tot_conc_label = kwargs.get('field_var', f'TOT_{dist_type.upper()}_CONC')
    tot_conc_title = kwargs.get('field_title', None)

    fig = plt.figure(figsize=plt.figaspect(.5))
    ax = fig.add_subplot(1, 2, 1)

    colors = dist_cmap(np.linspace(dist_cmap_normrange[0], 
                                   dist_cmap_normrange[1], 
                                   times.size))
    for c, time in zip(colors, times):
        x_vals = []
        bin_vals = []

        bin_edges = scenario_distdata['BIN_EDGES'][:].data[0]
        bin_width = bin_edges[1:] - bin_edges[:-1]

        for bin_idx in range(100):
            bin_idx += 1 # 1 indexing 
            bin_data = scenario_distdata[f'{dist_type}_a{str(bin_idx).zfill(3)}'][time, k, j, i].data.item()#/1e6
            if local_binning:
                bin_data = scenario_distdata[f'{dist_type}_a{str(bin_idx).zfill(3)}'][time, k, j-local_binning:j+local_binning, i-local_binning:i+local_binning].data#/1e6
                bin_data = bin_data.mean()
            bin_vals.append(bin_data)
            x_vals.append(bin_idx)

        ax.plot(scenario_distdata['BIN_CENTERS'][:].data[0], scaling_factor*np.array(bin_vals), 
                label=f'h {(1/60)*Archive.historydelta_m*time:1.0f}', c=c, lw=1.5)
        ax.set_xscale('log')

    ax.set_xlim(xlims[0], xlims[1])
    ax.set_ylim(ylims[0], ylims[1])
    ylims = ax.get_ylim()
    if yscale == 'log':
        ax.set_yscale('log')
    ax.legend(loc=legend_loc, handlelength=1, ncol=legend_ncol, columnspacing=0.8)

    if dist_type == 'num':
        ax.set_ylabel('Number concentration (m$^{-3}$)', fontsize=label_size)
        ax.set_title(f'Number distribution', fontsize=title_size)
    if dist_type == 'mass':
        ax.set_ylabel('Mass concentration ($\mu$g m$^{-3}$)', fontsize=label_size)
        ax.set_title(f'Mass distribution', fontsize=title_size)

    ax.set_xlabel('Diameter (m)', fontsize=label_size)

    ax.grid(which = "major", linewidth = 1, axis='y', ls="dashed", dashes=(4,4), c='#414141', alpha=.5)
    ax.grid(which = "minor", linewidth = 1, axis='y', ls="dashed", dashes=(6,6), c='white')
    ax.grid(which = "minor", linewidth = 1, axis='x', ls="dotted", dashes=(.5,6), c='#414141')
    ax.grid(which = "major", linewidth = 1, axis='x', ls="dotted", dashes=(.5,6), c='#414141')
    ax.tick_params(axis='both', labelsize=13, which='major', direction='in', top=True, right=True, bottom=True, left=True)
    ax.tick_params(axis='both', which='minor',direction='in',top=True, right=True, bottom=True, left=True)

    # Second subplot
    
    ax = fig.add_subplot(1, 2, 2, projection='3d')
    X, Y = np.meshgrid(np.arange(Archive.domain_x_cells), np.arange(Archive.domain_y_cells))

    
    conc_min = scenario_aerodata[tot_conc_label][totconctimeidx, 0, :, :].min()
    conc_max = scenario_aerodata[tot_conc_label][totconctimeidx, 0, :, :].max()

    if lognorm:
        norm = mplcolors.LogNorm(conc_min, conc_max)
    else:
        norm = None

    ax.contourf(X, Y, scenario_aerodata[tot_conc_label][totconctimeidx, k, :, :], 
                cmap=plt.cm.Spectral_r, norm=norm, zdir='z', offset=k)
    if local_binning:

        # plot a semi-transparent rectangle around the region of local binning
        data = scenario_aerodata[tot_conc_label][totconctimeidx, k, :, :].data
        data[:] = np.nan
        data[j-local_binning:j+local_binning, i-local_binning:i+local_binning] = .1
        
        ax.contourf(X, Y, data, 
                cmap=plt.cm.Greys, zdir='z', vmin=0, vmax=1, offset=k, alpha=.5)

        # plot a thin black boarder around the semi-transparent region
        width = 2
        data[j-local_binning-width:j-local_binning, i-local_binning-width:i+local_binning+width] = 1
        data[j+local_binning:j+local_binning+width, i-local_binning-width:i+local_binning+width] = 1
        data[j-local_binning:j+local_binning, i-local_binning-width:i-local_binning] = 1
        data[j-local_binning:j+local_binning, i+local_binning:i+local_binning+width] = 1
        data[j-local_binning:j+local_binning, i-local_binning:i+local_binning] = np.nan
        ax.contourf(X, Y, data, 
                cmap=plt.cm.Greys, zdir='z', vmin=0, vmax=1, offset=k, alpha=.8)
    else:
        ax.plot([i], [j], zs=k, zdir='z', marker=kwargs.get('subset_marker', '*'), c='white', markeredgecolor='k', 
            markeredgewidth=.5, markersize=12, zorder=10)

    nx = scenario_aerodata.dimensions['west_east'].size
    ny = scenario_aerodata.dimensions['south_north'].size
    nz = scenario_aerodata.dimensions['bottom_top'].size
    ax.set_xlim(0, nx)
    ax.set_ylim(0, ny)
    ax.set_zlim(0, nz)
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')

    xwidth = (nx*100) / 1000 # width of domain in km (assume standard cell width of 100 m)
    ywidth = (ny*100) / 1000 # width of domain in km (assume standard cell width of 100 m)
    zheight = 2 # height of domain in km (assume standard cell width of 100 m)

    ax.set_xticks(np.linspace(0, nx+1, 5))
    ax.set_xticklabels(np.linspace(0, xwidth, 5))
    ax.set_yticks(np.linspace(0, ny+1, 5))
    ax.set_yticklabels(np.linspace(0, ywidth, 5))
    ax.set_zticks(np.linspace(0, nz+1, 5))
    ax.set_zticklabels(np.linspace(0, zheight, 5))

    if tot_conc_label == f'TOT_{dist_type.upper()}_CONC':
        if dist_type == 'num':
            ax.set_title(f'Total number conc. ($t = ${Archive.historydelta_m*totconctimeidx/60} h)', y=1.17, 
                        transform=ax.transAxes, fontsize=title_size)
        if dist_type == 'mass':
            ax.set_title(f'Total mass conc. ($t = ${Archive.historydelta_m*totconctimeidx/60} h)', y=1.17, 
                        transform=ax.transAxes, fontsize=title_size)
    else:
        ax.set_title(f'{tot_conc_title} ($t = ${Archive.historydelta_m*totconctimeidx/60} h)', y=1.17, 
                        transform=ax.transAxes, fontsize=title_size)

    N = 1.5  # some number > 1 that stretches z axis as you desire
    ax.set_box_aspect((1, 1, N))  # xy aspect ratio is 1:1, but stretches z axis
    ax.view_init(elev=40., azim=-50)
    norm = matplotlib.colors.Normalize(vmin=conc_min, vmax=conc_max)
    sc = matplotlib.cm.ScalarMappable(cmap=plt.cm.Spectral_r, norm=norm)
    sc.set_array([])
    #cax = ax.inset_axes([.1,-.25, 0.8, 0.05]) #horizontal

    cax = ax.inset_axes([1.15,.0, 0.05, 0.9])
    plt.colorbar(sc, cax=cax, orientation='vertical', shrink=.75)

    if savefig:
        plt.savefig(f'{scenario}_{dist_type}conc_i{i}_j{j}_k{k}.pdf', format='pdf', bbox_inches='tight')
    plt.show()

    return bin_vals 
"""
def plotCCNError(variable, absolute=False):
    if variable not in boxplot_data:
        print(f'Calculating boxplot data for {variable}')
        calculateBoxplotData(variable)

    boxplot_vardata = boxplot_data[variable]
    if absolute:
        for i, data in enumerate(boxplot_vardata):
            boxplot_vardata[i] = abs(data)
    
    scenario_nsh_values = []
    boxplot_positions = []
    for scenario_name in emissions_nsh_dict.keys():
        scenario_nsh = emissions_nsh_dict[scenario_name]
        scenario_nsh_values.append(scenario_nsh)
        boxplot_positions.append(round(scenario_nsh, 2))

    # Creating boxplots
    medianprops = dict(linestyle='-', linewidth=1, color='#1B2860')  # Customize mean line color
    meanprops= dict(marker="o", markerfacecolor="#E9562B", markeredgecolor="#1B2860", markeredgewidth=.7)
    boxprops = dict(facecolor='#408BC5', edgecolor='#408BC5')
    whiskerprops = dict(color='#1B2860')
    width = 0.25
    fig, ax = plt.subplots(1,1, figsize=(7, 4))
    bp = ax.boxplot(boxplot_vardata, #positions=boxplot_positions, 
                    vert=True, patch_artist=True,
                    medianprops=medianprops, meanprops=meanprops, boxprops=boxprops,
                    whiskerprops=whiskerprops,
                    widths=width, showmeans=True,
                    showfliers=False) # dont display outliers

    # Adding labels and title
    x_labels = [f'{value:3.2f}' for value in scenario_nsh_values]
    #ax.set_xticks(group_centers)  # Positioning the ticks at the center of each group
    ax.set_xticklabels(x_labels)
    ax.set_xlabel('Emissions Spatial Heterogeneity', fontsize=13)
    if absolute:
        ylabel = r'Absolute % Error'
        absprefix = 'abs'
    else:
        ylabel = r'% Error'
        absprefix = ''
        ax.axhline(y=0, xmin=0, xmax=1, ls='--', c='k')
    ax.set_ylabel(ylabel, fontsize=13)
    if variable.startswith('ccn_'):
        variable_elements = variable.split('_')
        supersat =  int(variable_elements[1])/10
        title = f'CCN ($S={supersat:2.1f}\%$)'
    else:
        title = variable
    ax.set_title(title, fontsize=13)
    ax.tick_params(axis='both', labelsize=12)
    ax.invert_xaxis()

    # Add inset plots for emission patterns
    xtickslocs = ax.get_xticks()
    xlims = ax.get_xlim()
    norm_xticklocs = (xtickslocs-xlims[0])/(xlims[1]-xlims[0])
    subplot_dim = 0.16 # width and height of checkerboard subplot

    cwd = os.getcwd()
    shdir = 'sh-patterns'
    domain_x_cells = 40
    domain_y_cells = 40
    griddir = f'xres{domain_x_cells}yres{domain_y_cells}'
    basecase_filename = 'uniform-basecase.csv'
    basecase_array_path = os.path.join(cwd, shdir, griddir, basecase_filename)
    basecase_arr = np.genfromtxt(basecase_array_path, delimiter=',')
    cmap = plt.cm.get_cmap('Greys_r')
    for scenario_name, xtickloc in zip(emissions_nsh_dict.keys(), norm_xticklocs):
        
        filename = f'{scenario_name}.csv'
        array_path = os.path.join(cwd, shdir, griddir, filename)
        scenario_arr = np.genfromtxt(array_path, delimiter=',')
        scaling_factor = basecase_arr.sum() / scenario_arr.sum()
        scenario_arr = scaling_factor*scenario_arr
        
        # see https://stackoverflow.com/questions/42932908/matplotlib-coordinates-tranformation
        # for more on insetting subplots
        axins2 = il.inset_axes(ax, "100%", "100%", loc=3, borderpad=0,
            bbox_to_anchor=(xtickloc-0.5*subplot_dim,-0.4, subplot_dim, subplot_dim), 
            bbox_transform=ax.transAxes)
        vmin=1
        vmax=1e3
        axins2.pcolormesh(scenario_arr,#norm=mplcolors.LogNorm(vmin, vmax), 
                        cmap='Greys', vmin=-.35, vmax=1.5
                        )
        axins2.set_axis_off()
        axins2.set_aspect('equal')

    plt.savefig(f'{variable}-{absprefix}error-vs-nsh.pdf', format='pdf', bbox_inches='tight')
"""
