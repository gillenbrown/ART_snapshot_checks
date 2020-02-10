#include <stdio.h>
#include <string.h>
#include <libgen.h>
#include <math.h>
#include <inttypes.h>
#include "read_tree.h"

int merger(struct halo this_halo){
    // A satellite is accreted when the progenitor of a satellite's 
    // descendent is not itself.
    if (this_halo.desc != NULL && this_halo.prog != NULL){
        if (this_halo.desc->prog->id != this_halo.id){
            return 1;
        }
    }
    return 0;
}

int merger_in_main_tree(struct halo candidate, struct halo endpoint){
    // Check whether this merger merged directly with one of the progenitors
    // of the main halo, rather than some other halo that later merged with
    // the one we're interested in.
    // First check that it ended up in the halo of interest
    if (candidate.tree_root_id != endpoint.id){
        return 0;
    }
    // Then check that the descendant of the candidate is one of the 
    // progenitors of the endpoint halo
    struct halo endpoint_progenitor = endpoint;
    struct halo merger_result = *candidate.desc;
    while (1){
        if (endpoint_progenitor.id == merger_result.id){
            return 1; // found a match
        }
        else if (endpoint_progenitor.prog == NULL){
            return 0; // reached end of progenitors with no match
        }
        // advance our counter
        endpoint_progenitor = *endpoint_progenitor.prog;
    }
}

void write_mergers(struct halo satellite, FILE *out_file){
    // Here satellite should be a halo in it's last timestep of existence 
    // before being eaten by a larger halo.
    // The most representative mass of the mergers are the masses of the halos
    // as the satellite is falling in. Therefore we need to trace back to find
    // the most massive progenitor of our satellite. I still want to write when
    // it finally merged for completeness
    struct halo central = *satellite.desc->prog;
    struct halo satellite_progenitor = satellite;
    struct halo central_progenitor = central;
    float central_mass = 0;
    float satellite_mass = 0;
    float scale_factor = 0;
    int64_t central_id = 0;
    int64_t satellite_id = 0;
    while (1){
        if (satellite_progenitor.mvir > satellite_mass){
            central_mass = central_progenitor.mvir;
            satellite_mass = satellite_progenitor.mvir;
            scale_factor = satellite_progenitor.scale;
            central_id = central_progenitor.id;
            satellite_id = satellite_progenitor.id;
        }
        // check if we can go back a timestep
        if (satellite_progenitor.prog != NULL && 
            central_progenitor.prog != NULL){
            // If we can go back, put both of these equal to their 
            // next progenitor
            satellite_progenitor = *satellite_progenitor.prog;
            central_progenitor = *central_progenitor.prog;
        }
        else{ // We can't trace one of these halos back any further
            break;
        }
    }
    // Now we have the information we can write to the file.
    fprintf(out_file, "%15f %25E %25E %15lu %15lu %15f %25E %25E %15lu %15lu\n", 
            scale_factor, satellite_mass, central_mass, 
            satellite_id, central_id, central.scale, satellite.mvir,
	    central.mvir, satellite.id, central.id);
}

void all_accretions(struct halo endpoint, FILE *out_file){
    struct halo this_halo;
    for (int i=0; i < all_halos.num_halos; i++){
        this_halo = all_halos.halos[i];
        if (merger(this_halo)){
            if (merger_in_main_tree(this_halo, endpoint)){
                write_mergers(this_halo, out_file);
            }
        }
    }
}

void growth_history(struct halo endpoint, FILE *out_file){
    // start at the endpoint, and write out the progenitor at all redshifts
    struct halo endpoint_progenitor = endpoint;
    while (1){
        fprintf(out_file, "%15f %25E %15lu\n",
                endpoint_progenitor.scale,
                endpoint_progenitor.mvir,
                endpoint_progenitor.id);
        // if there is no more progenitor, end
        if (endpoint_progenitor.prog == NULL){
            return;
        }
        // otherwise, go to the next progenitor
        endpoint_progenitor = *endpoint_progenitor.prog;
    }
}

void write_satellite(struct halo central, struct halo satellite, 
                     FILE *out_file){
    // here we want to write the distance and peak mass of this satellite
    // distance is pretty easy. Note that this distance is in Mpccm/h
    float distance;
    distance = sqrt(pow(central.pos[0] - satellite.pos[0], 2) +
                    pow(central.pos[1] - satellite.pos[1], 2) +
                    pow(central.pos[2] - satellite.pos[2], 2));
    // turn distance into kpccm/h
    distance *= 1000;

    // we have to track the mass backwards to do get the peak mass
    struct halo satellite_progenitor = satellite;
    float satellite_mass = 0;
    while (1){
        if (satellite_progenitor.mvir > satellite_mass){
            satellite_mass = satellite_progenitor.mvir;
        }
        // check if we can go back a timestep
        if (satellite_progenitor.prog != NULL){
            // If we can go back, put it equal to their next progenitor
            satellite_progenitor = *satellite_progenitor.prog;
        }
        else{ // We can't trace one this back any further
            break;
        }
    }
    // Now we have the information we can write to the file.
    fprintf(out_file, "%30E %40f\n", 
            satellite_mass, distance);
}

void satellites(struct halo central, FILE *out_file){
    // first find all satellites of the central
    struct halo this_halo;
    for (int i=0; i < all_halos.num_halos; i++){
        this_halo = all_halos.halos[i];
        // check that they're at the last output
        if (this_halo.scale == central.scale){
            // the parent attribute tells what is a subhalo of what
            if (this_halo.parent != NULL){
                if (this_halo.parent->id == central.id){
                    write_satellite(central, this_halo, out_file);
                }
            }
        }
    }
}

char *make_output_files(char *tree_location, char *file_name){
    // Take the directory of the tree outputs and get the checks directory
    // inside the main directory
    // dirname returns the parent directory each time, so we can use it to
    // go up the number of levels to get to the main directory
    char temp_tree_loc[500];
    strcpy(temp_tree_loc, tree_location);
    static char temp_dir[500];
    strcpy(temp_dir, "");
    strcat(temp_dir, dirname(dirname(dirname(temp_tree_loc))));
    // then append the checks directory
    strcat(temp_dir, "/checks/");
    strcat(temp_dir, file_name);
    return temp_dir;
}

int main(int argc, char *argv[]) {

    // First we want to organize the output files, which we first need to do
    // Organize the output files
    FILE *mergers_1, *mergers_2, *growth_1, *growth_2, *satellites_1, *satellites_2;
    mergers_1 = fopen(make_output_files(argv[1], "mergers_1.txt"), "w");
    mergers_2 = fopen(make_output_files(argv[1], "mergers_2.txt"), "w");
    growth_1 = fopen(make_output_files(argv[1], "growth_1.txt"), "w");
    growth_2 = fopen(make_output_files(argv[1], "growth_2.txt"), "w");
    satellites_1 = fopen(make_output_files(argv[1], "satellites_1.txt"), "w");
    satellites_2 = fopen(make_output_files(argv[1], "satellites_2.txt"), "w");

    // First write the headers
    fprintf(mergers_1, "# %13s %25s %25s %15s %15s %15s %25s %25s %15s %15s\n",                            
            "scale_factor", "satellite_mass_[Msun/h]", "central_mass_[Msun/h]",   "satellite_id", "central_id",
            "final_scale",  "final_sat_mass_[Msun/h]", "final_cen_mass_[Msun/h]", "final_sat_id", "final_cen_id");
    fprintf(mergers_2, "# %13s %25s %25s %15s %15s %15s %25s %25s %15s %15s\n",                            
            "scale_factor", "satellite_mass_[Msun/h]", "central_mass_[Msun/h]",   "satellite_id", "central_id",
            "final_scale",  "final_sat_mass_[Msun/h]", "final_cen_mass_[Msun/h]", "final_sat_id", "final_cen_id");
    fprintf(growth_1, "# %13s %25s %15s\n", 
            "scale_factor", "central_mass_[Msun/h]", "central_id");
    fprintf(growth_2, "# %13s %25s %15s\n", 
            "scale_factor", "central_mass_[Msun/h]", "central_id");
    fprintf(satellites_1, "# %28s %40s\n", 
            "satellite_m_peak[Msun/h]", "satellite_final_distance_[kpccm/h]");
    fprintf(satellites_2, "# %28s %40s\n", 
            "satellite_m_peak[Msun/h]", "satellite_final_distance_[kpccm/h]");

    read_tree(argv[1]);

    // determine the scale factor of the last output
    float scale_last = 0;
    float this_scale;
    for (int i=0; i < all_halos.num_halos; i++){
        this_scale = all_halos.halos[i].scale;
        if (this_scale > scale_last){
            scale_last = this_scale;
	}
    }
    float comparison_scale = scale_last - 0.000001;
    printf("Last output: a = %f\nComp output: a = %f\n", scale_last, comparison_scale);

    // Get the biggest and smallest halos at the last output
    struct halo empty = {0};
    empty.mvir = 0;
    struct halo halo_1 = empty;
    struct halo halo_2 = empty;
    struct halo halo_n = empty;
    halo_n.mvir = 1E20;
    struct halo this_halo;
    for (int i=0; i < all_halos.num_halos; i++){
        this_halo = all_halos.halos[i];
        if (this_halo.scale > comparison_scale){
            if (this_halo.mvir > halo_1.mvir){
                halo_2 = halo_1;
                halo_1 = this_halo;
            }
            else if (this_halo.mvir > halo_2.mvir){
                halo_2 = this_halo;
            }
            else if (this_halo.mvir < halo_n.mvir){
                halo_n = this_halo;
            }
        }
    }
    // printf("Biggest Mass: %E\n", biggest.mvir);
    // printf("Second Mass: %E\n", second_biggest.mvir);
    printf("Smallest Mass: %E\n", halo_n.mvir);
    // write the accretion history files
    all_accretions(halo_1, mergers_1);
    all_accretions(halo_2, mergers_2);

    // then the growth history files
    growth_history(halo_1, growth_1);
    growth_history(halo_2, growth_2);

    // then write the satellite info
    satellites(halo_1, satellites_1);
    satellites(halo_2, satellites_2);

    fclose(mergers_1);
    fclose(mergers_2);
    fclose(growth_1);
    fclose(growth_2);
    fclose(satellites_1);
    fclose(satellites_2);

    return 0;
}
