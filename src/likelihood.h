#include "moma_input.h"


void sc_likelihood(const std::vector<double> &params_vec, MOMAdata &cell, double &total_likelihood){
    double likelihood = 0;
    if (cell.is_root()){
        /* use initial conditions for nm, nC */
        ;
    }
    else{
        /* update nm and nC that depend on mother cell */
        ;
    }

    for (int i=0; i<params_vec.size();++i){
        likelihood += (i+1)*params_vec[i];
    } 
    likelihood = pow(likelihood,2) + 1;

    cell.likelihood = likelihood;
    total_likelihood += likelihood;
    /* save likelihood, nm, nC in cell */
}


void likelihood_recr(const std::vector<double> &params_vec, MOMAdata *cell, double &total_likelihood){
    /*  
    * Recursive implementation that applies the function func to every cell in the genealogy
    * not meant to be called directly, see wrapper below
    */
    if (cell == nullptr)
        return;
    sc_likelihood(params_vec, *cell, total_likelihood);

    likelihood_recr(params_vec, cell->daughter1, total_likelihood);
    likelihood_recr(params_vec, cell->daughter2, total_likelihood);
}


double total_likelihood(const std::vector<double> &params_vec, std::vector<double> &grad, void *c){
    /*
    * total_likelihood of cell tree, to be maximized
    */

    // MOMAdata cell = *(MOMAdata *) c;
    double total_likelihood = 0;

    likelihood_recr(params_vec,  (MOMAdata *) c, total_likelihood);
    return total_likelihood;
}



/* ==========================================================================
* DEMO on how this work
*/

// example function to illustrate how this works
void set_generation(const std::vector<double> &params_vec, MOMAdata &cell){
    if (cell.parent != nullptr){
        cell.generation = cell.parent->generation + 1;
    } else{
        cell.generation = 0;
    }
}


void print_generation_tree(const std::vector<double> &params_vec, MOMAdata &cell){

    // apply the set_generation function to root cell followed by the first generation etc...
    apply_down_tree(params_vec, cell, set_generation);                    

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
