import numpy as np
import matplotlib.pyplot as plt

try:
    import cmasher as cm
    has_cmasher =  True
except ImportError:
    has_cmasher = False
    

plt.ion()

def plot_slice(array, vec_dim, slice_dim, shp, rfp, inc, vmin, vmax, show_cbar=True, show_labels=True, save_fig=False):
    """
    produces a plot showing a slice through GMF model

    :param array: GMF field in grid form e.g. from on_grid()
    :param vec_dim: which dimenion of vec(B) to display, can be 0, 1, 2. If set to any other value, the amplitude of B is taken.
    :param slice_dim: dimension on which to slice, e.g. 0 will give a slice through the y-z plane
    :param shp: shape of the image to plot, given as a list containing number of cells in [x, y, z]
    :param rfp: reference point = location of the smallest coordinate value in each direction, as list [refx, refy, refz]
    :param inc: increment in each dimension = size of cells used for plot, as list [incx, incy, incz]
    :param vmin: minimum value on color bar (magnetic field in G)
    :param vmax: maximum value on color bar (magnetic field in G)
    :param show_cbar: show color bar
    :param show_labels: show labels on plotted axes
    :param save_fig: saves figure automatically
    """
   
    if has_cmasher:
        cmap = getattr(cm, "prinsenvlag") 
    else:
        cmap = "RdBu_r"
    fig, ax = plt.subplots()
    
    slices = [slice(0, shp[0], None), slice(0, shp[1], None), slice(0, shp[2], None)]
    try:
        cut_index = int(shp[slice_dim]/2)
    except IndexError:
        raise IndexError('slice_dim must be 0, 1, or 2!')
    slices[slice_dim] = slice(cut_index, cut_index +1, 1)
    slices = tuple(slices)

    dims = [0, 1, 2]
    dims.remove(slice_dim)
    dims_label = ['x', 'y', 'z']
    try:
        comp_label = dims_label[vec_dim] + '_component'
        ax.imshow(np.squeeze(array[vec_dim][slices]).T, cmap=cmap, origin='lower', vmin=vmin, vmax=vmax)
    except IndexError:
        comp_label = 'amplitude'
        ax.imshow(np.squeeze(np.linalg.norm(array, axis=0)[slices]).T, cmap='Reds', origin='lower', vmin=0, vmax=vmax)
        
    dims_label.remove(dims_label[slice_dim])
    

    xticks_label = [int(rfp[dims[0]] + i*shp[dims[0]]/4*inc[dims[0]]) for i in range(5)]
    xticks_loc = [i*shp[dims[0]]/4 for i in range(5)]
    yticks_label = [int(rfp[dims[1]] + i*shp[dims[1]]/4*inc[dims[1]]) for i in range(5)]
    yticks_loc = [i*shp[dims[1]]/4 for i in range(5)]
        

    if show_labels:
        ax.set_yticks(yticks_loc, labels=yticks_label)
        ax.set_ylabel(r"$%s$ / kpc" % dims_label[1])
        ax.set_xticks(xticks_loc, labels=xticks_label)
        ax.set_xlabel(r"$%s$ / kpc" % dims_label[0])    
    else:
        plt.axis('off')


    if show_cbar:
        cbar = fig.colorbar(ax.images[0], orientation="horizontal", shrink=0.5, aspect=30)
        if show_labels:
            if comp_label == 'amplitude':
                cbar.set_label(r'$B$ / $\mu$G')
            else:
                cbar.set_label(r'$B_%s$ / $\mu$G' % comp_label[0])
    plt.tight_layout()
    if save_fig:
        plt.savefig(comp_label + "_" + dims_label[0] + "_" + dims_label[1] + "_plane")    
    plt.show()
    plt.close()