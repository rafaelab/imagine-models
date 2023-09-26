import numpy as np
import matplotlib.pyplot as plt

try:
    import cmasher as cm
    has_cmasher =  True
except ImportError:
    has_cmasher = False
    

plt.ion()

def plot_slice(array, vec_dim, slice_dim, shp, rfp, inc, vmin, vmax, show_cbar=True, show_labels=True, save_fig=False):
   
    if has_cmasher:
        cmap = getattr(cm, "prinsenvlag") 
    else:
        cmap = "RdBu_r" 
    fig, ax = plt.subplots()
    dims = [0, 1, 2]
    dims.remove(slice_dim)
    

    slices = [slice(0, shp[0], None), slice(0, shp[1], None), slice(0, shp[2], None)]
    cut_index = int(shp[slice_dim]/2)
    slices[slice_dim] = slice(cut_index, cut_index +1, 1)
    slices = tuple(slices)
    
    ax.imshow(np.squeeze(array[vec_dim][slices]).T, cmap=cmap, origin='lower', vmin=vmin, vmax=vmax)

    dims_label = ['x', 'y', 'z']
    comp_label = dims_label[vec_dim]
    dims_label.remove(comp_label)

    xticks_label = [int(rfp[dims[0]] + i*shp[dims[0]]/5*inc[dims[0]]) for i in range(6)]
    xticks_loc = [i*shp[dims[0]]/5 for i in range(6)]
    yticks_label = [int(rfp[dims[1]] + i*shp[dims[1]]/5*inc[dims[1]]) for i in range(6)]
    yticks_loc = [i*shp[dims[1]]/5 for i in range(6)]
        

    if show_labels:

        ax.set_yticks(yticks_loc, labels=yticks_label)
        ax.set_ylabel("kpc")
        ax.set_xticks(xticks_loc, labels=xticks_label)
        ax.set_xlabel("kpc")    
        ax.set_title(comp_label + " component in the  " + dims_label[0] + "-" + dims_label[1] + " plane")
    else:
        plt.axis('off')


    if show_cbar:
        fig.colorbar(ax.images[0], orientation="horizontal")
    plt.tight_layout()
    if save_fig:
        plt.savefig(comp_label + "_component_" + dims_label[0] + "_" + dims_label[1] + "_plane")    
    plt.show()
    plt.close()