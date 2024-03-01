from load_griddedoutput_datasets import *
from scipy import stats

def printGriddedOutputInfo(data):
    for var in data.variables.keys():
        try:
            descrip = data[var].long_name
        except AttributeError:
            try:
                descrip = data[var].description
            except AttributeError:
                descrip = 'no description'
        print(f'{var:27s}  dimensions: {str(data[var][:].shape):12s} description: {descrip}  -  ')

def getParticleDiameters(data):
    particle_volume = (data['aero_particle_mass'][:].T/data['aero_density'][:].reshape(1, 20)).sum(axis=1)
    particle_diameter = np.cbrt(particle_volume*(6/np.pi))
    return particle_diameter

def getGasData(data, species_name):
    # get the mixing ratio column data for passed gas species name
    gas_species_index = gas_species.index(species_name)
    return data['gas_mixing_ratio'][gas_species_index, :]

def getVertGridCellPartIndices(data, k):
    # return the start and end indices of the particles in the requested kth vertical 
    # grid cell
    # TODO: with this method I'm getting one less particle in the topmost grid cell
    # when I check against data['n_parts'][-1].item()
    #k -= 1 # convert from one indexing to zero indexing 

    start_idx = data['part_start_index'][k]
    n_parts_in_cell = data['n_parts'][k].item()
    end_idx = start_idx + n_parts_in_cell
    return start_idx, end_idx

def getNumberDist(particle_diam_arr, particle_numconc_arr, n_grid_cells=1):
    
    # returns array of size equal to the number of particles in the passed array
    # with elements equal to the bin number in which the particle is sorted into
    digitized = np.digitize(particle_diam_arr, bin_edges) 
        
    # one could use the digitized array as a map to access particles within a certain 
    # bin number range (below I'm getting particles in bins below #40 but you could
    # do both lower and upper bounds). A potentially useful extension of this would be
    # to pass upper and lower bounds in terms of particle diameters, for which the closest
    # bin edges are found and corresponding indicies, which would then be passed on to this
    # mapping as lower/upper bounds on the digitized array
    #sfc_part_diams[digitized<40]

    # OLD METHOD OF CALCULATING NUMBER DIST, SEE CODE BELOW FOR UPDATED ROUTINE
    # get the number of particles in each bin
    #bin_count = np.bincount(digitized, minlength=bin_edges.size)
    #dNdlogDp = bin_count[:-1] / bin_logwidth 

    # My hunch is that my current way of plotting the number distribution is somewhat 
    # incorrect given that im not accounting for the computational particle multiplicity. 
    # I ought to be using the `aero_num_conc` variable, which indicates the number concentration 
    # for each computational particle in a grid cell. 
    # For each bin, find the total number concentration by using the digitized array as a mapping
    bin_num_conc = np.zeros(n_bins)
    for bin_idx in np.arange(n_bins):
        # note the need to subtract 1 from digitized since the bin values start at 1
        bin_num_conc[bin_idx] = particle_numconc_arr[(digitized-1)==bin_idx].sum()
    dNdlogDp = bin_num_conc[:-1] / bin_logwidth
    dNdlogDp = dNdlogDp/n_grid_cells

    #print(bin_count.sum()) # verify the total bin count matches number of particles

    return dNdlogDp

def getNumberDistOptimized(particle_diam_arr, particle_numconc_arr, n_grid_cells=1):
    """Used scipy's binned_statisic method to efficiently sort particles by diameter 
    and compute number concentration 

    Jeff Curtis recommended using this method instead of digitized since its a bit more
    compact and likely more efficient
    """
    values = particle_numconc_arr
    sorting_variable = particle_diam_arr
    # generate weighted binning of particles by diameter - sorted by diameter and the sum of number conc in 
    # each bin is the resulting histogram value
    bin_result, _ , binnumbers = stats.binned_statistic(sorting_variable, values, 'sum', bins=bin_edges)
    dNdlogDp = bin_result/bin_logwidth
    dNdlogDp = dNdlogDp/n_grid_cells
    return dNdlogDp

def getBinnedSpeciesMass(particle_diam_arr, particle_mass_arr, particle_numconc_arr, n_grid_cells=1):

    digitized = np.digitize(particle_diam_arr, bin_edges) 

    # create an array of size num species x number of bins
    # for each particle in particle diameter array
    # add amount of each species in particle to species in bin
    binned_species_mass = np.zeros((20, n_bins)) 
    for binidx, particle, num_conc in zip(digitized, particle_mass_arr.T, particle_numconc_arr):
        # note the need to subtract 1 from digitized binidx since the bin values start at 1
        binned_species_mass[:, binidx-1] += particle*num_conc
    
    dMdlogDp = binned_species_mass/bin_logwidth
    dMdlogDp = dMdlogDp/n_grid_cells # compute mean if particle array is across multiple grid cells

    return dMdlogDp

def getBinnedSpeciesMassOptimized(particle_diam_arr, particle_mass_arr, particle_numconc_arr, n_grid_cells=1):
    """Used scipy's binned_statisic method to efficiently sort particles by diameter 
    and compute mass concentration 

    Jeff Curtis recommended using this method instead of digitized since its a bit more
    compact and likely more efficient
    """

    sorting_variable = particle_diam_arr
    values = particle_mass_arr*particle_numconc_arr

    bin_result, _ , binnumbers = stats.binned_statistic(sorting_variable, values, 'sum', bins=bin_edges)
    dMdlogDp = bin_result/bin_logwidth
    dMdlogDp = dMdlogDp/n_grid_cells # compute mean if particle array is across multiple grid cells

    return dMdlogDp