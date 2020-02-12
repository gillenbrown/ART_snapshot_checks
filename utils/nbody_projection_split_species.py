import yt
import cmocean
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
import matplotlib.patheffects as PathEffects
import betterplotlib as bpl

bpl.presentation_style()

cmap = cmocean.cm.deep_r
# =========================================================================
#         
# Mega plot with all the individual N-body species
# 
# =========================================================================
# This is some ugly code here, but it was actually the easiest way to get a 
# multipanel plot. yt is not easy to handle, but this worked

def nbody_projection_split_species(ds, center, plot_width, axes_unit, savename):
    # Make a box so we don't have to project through the full simulation box
    left_edges = [c - (plot_width/2.0) for c in center]
    right_edges = [c + (plot_width/2.0) for c in center]
    box = ds.region(center=center, left_edge=left_edges, right_edge=right_edges)

    # then get the number of density fields and figure out how many panels 
    # we need the plot to have
    density_fields = [item for item in ds.derived_field_list
                      if item[0] == "deposit" 
                      and "N-BODY" in item[1] and "density" in item[1]]
    n_panels = len(density_fields)
    n_rows = min(2, n_panels)
    n_cols = int(np.ceil(n_panels / n_rows))
    extra_plot = (n_panels % n_rows) != 0  # get rid of the last one?

    # the labeling of this will be ugly, since we do things in terms of pixels
    width_tuple = (float(plot_width.to("kpc").value), "kpc")
    n_pix = 1E4
    # determine the mapping from pixels to Mpc
    pix_per_kpc = n_pix / width_tuple[0]
    center_pix = n_pix / 2.0

    # figure out the right value for the increments of the ticks on the axes
    base_kpc = 10**np.floor(np.log10(width_tuple[0]/2.0))
    base_pix = base_kpc * pix_per_kpc

    # start from the center and work our way outwards
    pix_vals = [center_pix]
    pix_labels = ["0"]
    i = 0
    while True:
        next_hi = center_pix + i*base_pix
        next_lo = center_pix - i*base_pix
        
        if next_hi < n_pix:
            pix_vals.append(next_hi)
            pix_vals.append(next_lo)

            pix_labels.append("{:g}".format(i*base_kpc))
            pix_labels.append("{:g}".format(-i*base_kpc))

            i += 1
        else:
            break

    # Create the projections. We can do all fields at once, then handle them later
    proj = yt.ProjectionPlot(ds, 'x', density_fields, width=plot_width,
                             center=center, data_source=box)
    # get the underlying data which we can use in imshow
    proj_data = proj.data_source.to_frb(width=width_tuple, height=width_tuple, 
                                        resolution=n_pix)

    # handle the colormaps
    norm = LogNorm(vmin=1E-6, vmax=0.1)
    # since we are logging the data, zeros are bad. Handle those in the colormap
    cmap.set_bad(cmap(0))

    # then actually make the plot
    fig = plt.figure(figsize=[8*n_cols, 8*n_rows])
    # have a grid of axes, plus one extra column for the colorbar
    gs = gridspec.GridSpec(nrows=n_rows, ncols=n_cols+1,
                           width_ratios=[1]*n_cols + [0.1], 
                           wspace=0.02, hspace=0.02)

    axs = []
    for r in range(n_rows):
        for c in range(n_cols):
            axs.append(fig.add_subplot(gs[r,c], projection="bpl"))
    # last column is the colorbar
    cax = fig.add_subplot(gs[:,n_cols], projection="bpl")

    # then go through each of these and plot the right data
    for i, field in enumerate(density_fields):  
        im_data = np.array(proj_data[field])
        im = axs[i].imshow(im_data, origin="lower", norm=norm, cmap=cmap)
        
        if field[1] == "N-BODY_density":
            title = "N-Body Total"
        else:
            title = field[1].replace("BODY", "Body")
            title = title.replace("_", " ")
            title = title.replace("density", "")
        
        text = axs[i].easy_add_text(title, "upper left", c="w")
        text.set_path_effects([PathEffects.withStroke(linewidth=5,
                               foreground=bpl.almost_black)])
        
    for i, ax in enumerate(axs):
        # set the limits, may be removed later
        ax.xaxis.set_ticks(pix_vals)
        ax.yaxis.set_ticks(pix_vals)
        ax.xaxis.set_ticklabels(pix_labels)
        ax.yaxis.set_ticklabels(pix_labels)    
        # remove the outer boxes
        ax.remove_spines(["all"])
        
        # then remove the plot labels depending on where we are
        row_number = i // n_cols
        col_number = i % n_cols
        
        # if we're not in the last row, we need to remove the x labels
        if row_number != (n_rows - 1):
            # but if we're the first item, we need to keep the y 
            if col_number == 0:
                ax.remove_labels("x")
                ax.add_labels(y_label="Y [kpc]")
            else:
                # otherwise we remove both
                ax.remove_labels("both")
        else:  # last row, keep x labels
            ax.add_labels(x_label="X [kpc]")   
            # remove y label unless first column
            if col_number == 0:
                ax.add_labels(y_label="Y [kpc]")
            else:
                ax.remove_labels("y")
        
    # make the colorbar
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label("Projected Density [g/cm$^2$]")

    # put the time and redshift at the top
    z = ds.current_redshift
    t = ds.current_time.to("Gyr").value
    fig.suptitle("t = {:.2f} Gyr, z = {:.2f}".format(t, z),
                 weight="bold", fontsize=24)

    # if we have an odd number of panels remove the last one
    if extra_plot:
        axs[-1].set_axis_off()

    fig.savefig(savename, dpi=400, bbox_inches='tight')
