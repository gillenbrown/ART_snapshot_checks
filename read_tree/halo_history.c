#include <stdio.h>
#include <string.h>
#include <libgen.h>
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
    // the most massive progenitor of our satellite.
    struct halo satellite_progenitor = *satellite.prog;
    struct halo central_progenitor = *satellite.desc->prog;
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
    fprintf(out_file, "%15f %15E %15E %15lu %15lu\n", 
            scale_factor, satellite_mass, central_mass, 
            satellite_id, central_id);
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
        fprintf(out_file, "%15f %15E %15lu\n",
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
    FILE *mergers_1, *mergers_2, *growth_1, *growth_2;
    mergers_1 = fopen(make_output_files(argv[1], "mergers_1.txt"), "w");
    mergers_2 = fopen(make_output_files(argv[1], "mergers_2.txt"), "w");
    growth_1 = fopen(make_output_files(argv[1], "growth_1.txt"), "w");
    growth_2 = fopen(make_output_files(argv[1], "growth_2.txt"), "w");

    // First write the headers
    fprintf(mergers_1, "# %13s %15s %15s %15s %15s\n", "scale_factor", 
            "satellite_mass", "central_mass", "satellite_id", "central_id");
    fprintf(mergers_2, "# %13s %15s %15s %15s %15s\n", "scale_factor", 
            "satellite_mass", "central_mass", "satellite_id", "central_id");
    fprintf(growth_1, "# %13s %15s %15s\n", 
            "scale_factor", "central_mass", "central_id");
    fprintf(growth_2, "# %13s %15s %15s\n", 
            "scale_factor", "central_mass", "central_id");

    read_tree(argv[1]);

    // Find all halos at z=0, then get the biggest ones
    struct halo empty = {0};
    empty.mvir = 0;
    struct halo halo_1 = empty;
    struct halo halo_2 = empty;
    struct halo this_halo;
    for (int i=0; i < all_halos.num_halos; i++){
        this_halo = all_halos.halos[i];
        if (this_halo.scale > 1.0){
            if (this_halo.mvir > halo_1.mvir){
                halo_2 = halo_1;
                halo_1 = this_halo;
            }
            else if (this_halo.mvir > halo_2.mvir){
                halo_2 = this_halo;
            }
        }
    }
    // printf("Biggest Mass: %E\n", biggest.mvir);
    // printf("Second Mass: %E\n", second_biggest.mvir);

    // write the accretion history files
    all_accretions(halo_1, mergers_1);
    all_accretions(halo_2, mergers_2);

    // then the growth history files
    growth_history(halo_1, growth_1);
    growth_history(halo_2, growth_2);

    fclose(mergers_1);
    fclose(mergers_2);
    fclose(growth_1);
    fclose(growth_2);

    return 0;
}
