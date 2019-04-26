#include <stdio.h>
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
    if (this_halo.tree_root_id != endpoint.id){
        return 0;
    }
    // Then check that the descendant of the candidate is one of the 
    // progenitors of the endpoint halo
    struct halo endpoint_progenitor = endpoint;
    struct halo merger_result = candidate.desc;
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
    while (1){
        if (satellite_progenitor.mvir > satellite_mass){
            central_mass = central_progenitor.mvir;
            satellite_mass = satellite_progenitor.mvir;
            scale_factor = satellite_progenitor.scale;
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
    // DO THIS NEXT
}

void all_accretions(struct halo endpoint){
    struct halo this_halo;
    for (int i=0; i < all_halos.num_halos; i++){
        this_halo = all_halos.halos[i];
        if (merger(this_halo)){
            if (merger_in_main_tree(this_halo, endpoint)){
                ;
            }
        }
    }
}

int main(int argc, char *argv[]) {
  printf(argv[1]);
  read_tree(argv[1]);
  printf("\n%"PRId64" halos found in tree_0_0_0.dat!\n", all_halos.num_halos);

  // Find all halos at z=0, then get the biggest ones
  struct halo empty = {0};
  empty.mvir = 0;
  struct halo biggest = empty;
  struct halo second_biggest = empty;
  struct halo this_halo;
  for (int i=0; i < all_halos.num_halos; i++){
    this_halo = all_halos.halos[i];
    if (this_halo.scale > 1.0){
        if (this_halo.mvir > biggest.mvir){
            second_biggest = biggest;
            biggest = this_halo;
        }
        else if (this_halo.mvir > second_biggest.mvir){
            second_biggest = this_halo;
        }
    }
  }
  printf("Biggest Mass: %E\n", biggest.mvir);
  printf("Second Mass: %E\n", second_biggest.mvir);

  // A satellite is accreted when the progenitor of a satellite's descendent is not itself.
  all_accretions(biggest);
  all_accretions(second_biggest);

  return 0;
}
