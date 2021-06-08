grid_level_field = ("index", "grid_level")


def level_projection(ds, center, plot_width, axes_unit, savename):
    grid_plot = yt.ProjectionPlot(
        ds, "x", grid_level_field, method="mip", center=center, width=plot_width
    )
    grid_plot.set_log(grid_level_field, False)
    grid_plot.set_cmap(grid_level_field, "tab20")
    grid_plot.set_zlim(grid_level_field, -0.5, 19.5)
    grid_plot.annotate_timestamp(
        redshift=True,
        corner="upper_left",
        time_unit="Gyr",
        time_format="t = {time:.2f} {units}",
        redshift_format="z = {redshift:.2f}",
    )
    grid_plot.set_axes_unit(axes_unit)

    grid_plot.save(savename, mpl_kwargs={"dpi": 400})
