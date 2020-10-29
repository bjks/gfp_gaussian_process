#include "moma_input.h"

void sc_likelihood(MOMAdata &cell, Parameter_set& params){
    double likelihood = 0;
    if (cell.is_root()){
        /* code */
        ;
    }
    else if (cell.is_leaf()){
        /* code */
        ;
    }
    else{
        ;
    }
    
    for(long t=0; t<cell.time.size(); ++t){
        ;
    }
    /* save likelihood, nm, nC in cell */
}

double total_likelihood(MOMAdata &cell, Parameter_set& params){
    /*
    * total_likelihood of cell tree, to be maximized
    */

    apply_down_tree(cell, sc_likelihood, params);

    double total_likelihood=0;
    /*
    * Add likelihoods
    */
   return total_likelihood;
}


/*
* DEMO on how this work
*/

// example function to illustrate how this works
void set_generation(MOMAdata &cell, Parameter_set &params){
    if (cell.parent != nullptr){
        cell.generation = cell.parent->generation + 1 * params.var_x;
    } else{
        cell.generation = 0;
    }
}


void print_generation_tree(MOMAdata &cell, Parameter_set& params){

    // apply the set_generation function to root cell followed by the first generation etc...
    apply_down_tree(cell, set_generation, params);


    /* output the genealogy with generation
        -> 20150630.5.4.124 generation: 0 -> 20150630.5.4.133 generation: 1 
        -> 20150630.5.4.140 generation: 2 -> 20150630.5.4.148 generation: 3
    */
    std::vector<std::vector<MOMAdata *> > cell_paths = get_genealogy_paths(cell);

    std::cout << std::endl;
    for (std::vector<MOMAdata *> path : cell_paths){
        for (MOMAdata * cell : path){
            std::cout << " -> " << cell->cell_id << " generation: " << cell->generation ;
        }
        std::cout << "\n\n";
    }

}
