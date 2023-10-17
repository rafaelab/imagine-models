import numpy as np
import matplotlib.pyplot as plt

try:
    import cmasher as cm
    has_cmasher =  True
except ImportError:
    has_cmasher = False
    

plt.ion()

def plot_slice(array, vec_dim, slice_dim, shp, rfp, inc, vmin, vmax, show_cbar=True, show_labels=True,
               save_fig=False, field_name=None, quiver=False, amplitude=False, cut_index =None):
    """
    produces a plot showing a slice through a IMAGINE model output. The slice goes throgh the middle of array in the respective dimension slice_dim

    :param array: vector or scalar field in grid form e.g. from on_grid()
    :param vec_dim: which dimension of vec(B) to display, can be 0, 1, 2. No effect for scalar field or if amplitude=True
    :param slice_dim: dimension on which to slice, e.g. 0 will give a slice through the y-z plane
    :param shp: shape of the image to plot, given as a list containing number of cells in [x, y, z]
    :param rfp: reference point = location of the smallest coordinate value in each direction, as list [refx, refy, refz]
    :param inc: increment in each dimension = size of cells used for plot, as list [incx, incy, incz]
    :param vmin: minimum value on color bar (magnetic field in muG)
    :param vmax: maximum value on color bar (magnetic field in muG)
    :param show_cbar: show color bar
    :param show_labels: show labels on plotted axes
    :param save_fig: saves figure automatically
    :param field_name: if not None: displays name of field in plot
    :param quiver: adds field vectors displayed as arrows to plot, has no effect for scalar fields
    """
   #
    is_vector = False
    if isinstance(array, list): 
        is_vector = True 
   
    if amplitude:
        cmap = "Reds"
        vmin = 0 
    else:
        if has_cmasher:
            cmap = getattr(cm, "prinsenvlag") 
        else:
            cmap = "RdBu_r"
    
    fig, ax = plt.subplots()
    
    slices = [slice(0, shp[0], None), slice(0, shp[1], None), slice(0, shp[2], None)]
    if cut_index is None:
        try:
            cut_index = int(shp[slice_dim]/2)
        except IndexError:
            raise IndexError('slice_dim must be 0, 1, or 2!')
    slices[slice_dim] = slice(cut_index, cut_index +1, 1)
    slices = tuple(slices)

    dims = [0, 1, 2]
    dims.remove(slice_dim)
    dims_label = ['x', 'y', 'z']
    
    
    if is_vector:
        if amplitude:
            array = np.squeeze(np.linalg.norm(array, axis=0)[slices]).T
            comp_label = 'amplitude'
        else:
            comp_label = dims_label[vec_dim] + '_component'
            array = array[vec_dim]
        
    ax.imshow(np.squeeze(array[slices]).T, cmap=cmap, origin='lower', vmin=vmin, vmax=vmax)

    if is_vector:
        if quiver:  # add arrows indicating B-field direction
            coord1 = [i for i in range(shp[dims[0]])][::5]  # only every fifth cell gets an arrow for better overview
            coord2 = [i for i in range(shp[dims[1]])][::5]
            B1 = np.squeeze(array[dims[0]][slices])[::5, ::5]
            B2 = np.squeeze(array[dims[1]][slices])[::5, ::5]
            plt.quiver(coord1, coord2, B1.T, B2.T, pivot='mid')
    
    slice_dim_label = dims_label[slice_dim]
    dims_label.remove(slice_dim_label)
    
    if rfp[dims[0]] / float(int(rfp[dims[0]])) - 1 < 1e-3:
        xticks_label = [int(rfp[dims[0]] + i*shp[dims[0]]/4*inc[dims[0]]) for i in range(5)]
    else:
        xticks_label = ['%.2f' % (rfp[dims[0]] + i*shp[dims[0]]/4*inc[dims[0]]) for i in range(5)]
    xticks_loc = [i*shp[dims[0]]/4 for i in range(5)]

    if rfp[dims[1]] / float(int(rfp[dims[1]])) - 1 < 1e-3:
        yticks_label = [int(rfp[dims[1]] + i*shp[dims[1]]/4*inc[dims[1]]) for i in range(5)]
    else:
        yticks_label = ['%.2f' % (rfp[dims[1]] + i*shp[dims[1]]/4*inc[dims[1]]) for i in range(5)]
    yticks_loc = [i*shp[dims[1]]/4 for i in range(5)]
        

    if show_labels:
        ax.set_yticks(yticks_loc, labels=yticks_label)
        ax.set_ylabel(r"$%s$ / kpc" % dims_label[1])
        ax.set_xticks(xticks_loc, labels=xticks_label)
        ax.set_xlabel(r"$%s$ / kpc" % dims_label[0])
        keyword = '%s = %.3f kpc' % (slice_dim_label, rfp[slice_dim]+cut_index*inc[slice_dim]) 
        keyword = field_name + ', ' + keyword if field_name is not None else keyword
        ax.text(0, 1, keyword, va='bottom', transform = ax.transAxes)
    else:
        plt.axis('off')


    if show_cbar:
        cbar = fig.colorbar(ax.images[0], orientation="horizontal", shrink=0.5, aspect=30)
        if show_labels:
            if is_vector:
                if amplitude:
                    cbar.set_label(r'$B$ / $\mu$G')
                else:
                    cbar.set_label(r'$B_%s$ / $\mu$G' % comp_label[0])
            else: 
                cbar.set_label(r'$n_\mathrm{th}$ / $\mathrm{cm}^-3$')
    plt.tight_layout()
    if save_fig:
        name = comp_label + "_" + dims_label[0] + "_" + dims_label[1] + "_plane"
        if not is_vector:
            name = name[1:]
        plt.savefig(name)    
    plt.show()
    plt.close()