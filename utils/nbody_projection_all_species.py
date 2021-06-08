import cmocean
import yt

cmap = cmocean.cm.deep_r


def nbody_projection_all_species(ds, center, plot_width, axes_unit, savename):
    left_edges = [c - (plot_width / 2.0) for c in center]
    right_edges = [c + (plot_width / 2.0) for c in center]
    box = ds.region(center=center, left_edge=left_edges, right_edge=right_edges)

    n_body_density_field = ("deposit", "N-BODY_density")

    plot = yt.ProjectionPlot(
        ds,
        "x",
        n_body_density_field,
        center=center,
        width=plot_width,
        data_source=box,
        axes_unit=axes_unit,
    )
    plot.annotate_timestamp(
        redshift=True,
        corner="upper_left",
        time_unit="Gyr",
        time_format="t = {time:.2f} {units}",
        redshift_format="z = {redshift:.2f}",
    )

    plot.set_background_color(n_body_density_field, cmap(0))
    plot.set_cmap(n_body_density_field, cmap)
    plot.set_zlim(n_body_density_field, 1e-5, 0.1)

    if savename is not None:
        plot.save(savename, mpl_kwargs={"dpi": 400})
    else:
        return plot
