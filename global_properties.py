import sys
import os

import yt
import numpy as np
import matplotlib

import betterplotlib as bpl
yt.funcs.mylog.setLevel(50)  # ignore yt's output
bpl.presentation_style()
bpl.presentation_style()  # for some reason this needs to be there twice

def print_and_write(info, file_obj):
    file_obj.write(info + "\n")
    print(info)

if __name__ == "__main__":
    ds_loc = sys.argv[1]
    ds = yt.load(ds_loc)
    ad = ds.all_data()
    name = os.path.abspath(ds_loc).replace(os.sep, "_").replace(".art", "")
    name = name.lstrip("_")

    # get the location of where to write the file.
    file_dir = "/u/home/gillenb/analysis/snapshot_checks/"
    file_path = file_dir + "summary_" + name + ".txt"
    plots_dir = file_dir + "plots/"
    print("Output being written to:")
    print(file_path)

    # see if there is an existing file here that we don't want to replace.
    if os.path.isfile(file_path):
        good = False
        while not good:
            choice = input("This output file already exists. " +
                           "Overwrite? (y/n): ")
            if choice.lower() in ["y", "n"]:
                good = True
            else:
                print("Choose 'y' or 'n'")

        if choice.lower() == "n":
            # don't overwrite, so the code ends
            print("File won't be overwritten, so nothing will be done.")
            exit()
        # Don't need to do anything here inside this block if the user does 
        # want to overwrite, since the file created later will automatically 
        # overwrite any existing file.
    # open the file
    out_file = open(file_path, "w")

    # =========================================================================
    #         
    # Basic information
    # 
    # =========================================================================
    print()  # no write since this is just for console formatting
    print_and_write("Simulation location:", out_file)
    print_and_write(ds_loc, out_file)

    a = ds.scale_factor
    z = 1.0 / a - 1.0
    print_and_write("\na = {:.4f}".format(a), out_file)
    print_and_write("z = {:.4f}".format(z), out_file)

    # =========================================================================
    #         
    # Grid Structure
    # 
    # =========================================================================
    box_size = ds.domain_width.to("Mpccm")[0]
    full_grid_idxs = ad[('index', 'grid_level')].value
    grid_levels, num_in_grid = np.unique(full_grid_idxs, return_counts=True)
    total_cells = np.sum(num_in_grid)
    cell_sizes = np.unique(ad["index", "dx"]).to("kpccm")[::-1]
    # ^ np.unique returns the unique values in sorted order. We want the
    # largest cells to correspond to the smallest level, so we reverse it
    min_cell_size = np.min(cell_sizes)
    max_cell_size = np.max(cell_sizes)
    base_grid_size = np.log2(box_size / max_cell_size)

    print_and_write("\nGrid Structure\n=============", out_file)
    print_and_write("box size: {:.3f}".format(box_size), out_file)
    print_and_write("base grid size: {:.3f}".format(base_grid_size), out_file)
    print_and_write("total cells: {:,}".format(total_cells), out_file)
    grid_out_str = "{:<5d}    {:>10,}    {:<8.3f}    {:<8.3f}"
    grid_header_str = "{:<5}    {:>10}    {:<14}    {:<20}"
    print_and_write(grid_header_str.format("Level", "Num Cells", "Cell Size",
                                           "Expected Cell Size"), 
                    out_file)
    for level in grid_levels:
        level = int(level)
        num_cells = num_in_grid[level]
        cell_size = cell_sizes[level]
        expected_cell_size = box_size / (2**(level + base_grid_size))
        print_and_write(grid_out_str.format(level, num_cells, cell_size, 
                                            expected_cell_size.to("kpccm")),
                        out_file)

    # =========================================================================
    #         
    # Particles
    # 
    # =========================================================================
    all_masses = ad[('N-BODY', 'MASS')].to("Msun")
    masses, num_particles = np.unique(all_masses, return_counts=True)
    # we reverse for the same reason we reversed cell sizes
    masses = masses[::-1] 
    num_particles = num_particles[::-1]
    total_particles = np.sum(num_particles)

    print_and_write("\nParticles\n=========", out_file)
    print_and_write("total particles: {:,}".format(total_particles), out_file)
    part_out_str = "{:<5d}    {:>15,}    {:<.2e}    {:<.8f}"
    part_out_str_top = "{:<5d}    {:>15,}    {:<.2e}    ----------"
    part_header_str = "{:<5}    {:>15}    {:<13}    {:<25}"
    print_and_write(part_header_str.format("Level", "Num Particles", "Mass", 
                                           "Mass Ratio to Previous"), 
                    out_file)
    for level in range(len(masses)):
        level = int(level)
        dm_mass = masses[level]
        dm_num = num_particles[level]
        
        if level == 0:
            print_and_write(part_out_str_top.format(level, dm_num, dm_mass),
                            out_file)
        else:
            ratio = (masses[level - 1] / dm_mass).value
            print_and_write(part_out_str.format(level, dm_num, dm_mass, ratio),
                            out_file)

    # =========================================================================
    #         
    # Plots
    # 
    # =========================================================================
    grid_plot_name = plots_dir + "grid_idxs_{}.png".format(name)
    n_body_plot_name = plots_dir + "n_body_{}.png".format(name)
    print_and_write("\nPlots will be saved to:", out_file)
    print_and_write(grid_plot_name, out_file)
    print_and_write(n_body_plot_name, out_file)
    
    n_body_field = ("deposit", "N-BODY_density")
    grid_level_field = ('index', 'grid_level')
    
    grid_plot = yt.SlicePlot(ds, normal=[1, 0, 0], fields=grid_level_field, 
                             width=(15, "Mpccm"))
    grid_plot.set_log(grid_level_field, False)
    grid_plot.set_cmap(grid_level_field, "Pastel1")
    grid_plot.set_zlim(grid_level_field, -0.5, 8.5)
    grid_plot.save(grid_plot_name)

    n_body_plot = yt.SlicePlot(ds, normal=[1, 0, 0], fields=n_body_field, 
                               width=(15, "Mpccm"))
    n_body_plot.save(n_body_plot_name)

    out_file.close()