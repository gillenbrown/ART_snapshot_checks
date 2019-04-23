#include <stdio.h>
#include <inttypes.h>
#include "read_tree.h"

void all_accretions(struct halo endpoint){
    struct halo this_halo;
    for (int i=0; i < all_halos.num_halos; i++){
        this_halo = all_halos.halos[i];
        // Only pick the halos that are part of the given halo tree.
        if (this_halo.tree_root_id == endpoint.id){
            // A satellite is accreted when the progenitor of a satellite's descendent is not itself.
            if (this_halo.desc != NULL && this_halo.prog != NULL){
                if (this_halo.desc->prog->id != this_halo.id && this_halo.mvir > 1E8){
                  printf("%lu %f %E\n", this_halo.id, this_halo.scale, this_halo.mvir);
                }
                struct halo *desc = this_halo.desc;
                // struct halo *desc_parent = desc->parent;
                // printf("This:%lu descendent:%lu progenitor:%lu this scale:%f\n", this_halo.id, this_halo.desc->id, this_halo.prog->id, this_halo.scale);
            }
            
            // if (this_halo.desc->parent->id != this_halo.id){
            //     printf("%lu %f\n", this_halo.id, this_halo.scale);
            // }
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
