#include "read_config.h"
#include "moma_input.h"
#include "utils.h"
#include "likelihood.h"



#include <iostream> 
#include <iterator> 
#include <map> 


int main(int argc, char** argv){

    Config config("config.txt");
    std::string infile = argv[1];
    
    if(! std::__fs::filesystem::exists(infile)){
        std::cout << "File " << infile << " not found! \nQuit" << std::endl;
        return 0;
    }
    // Read data
    std::vector<MOMAdata> cells =  getData(infile, 
                                        config.time_col,
                                        config.length_col,
                                        config.fp_col,
                                        config.delm);
    build_cell_genealogy(cells);

    // get the "tree" starting from all root cells
    std::vector<MOMAdata *> root_cells = get_roots(cells);

    for(long j=0; j<root_cells.size(); ++j){
        std::vector<std::vector<MOMAdata *> > cell_paths = get_genealogy_paths(*root_cells[j]);

        // apply the set_generation function to root cell followed by the first generation etc...
        apply_down_tree(*root_cells[j], set_generation);


        /* output the genealogy with generation
         -> 20150630.5.4.124 generation: 0 -> 20150630.5.4.133 generation: 1 
            -> 20150630.5.4.140 generation: 2 -> 20150630.5.4.148 generation: 3
        */
        std::cout << std::endl;
        for (std::vector<MOMAdata *> path : cell_paths){
            for (MOMAdata * cell : path){
                std::cout << " -> " << cell->cell_id << " generation: " << cell->generation ;
            }
            std::cout << "\n\n";
        }

    }
    
    std::cout << "Done." << std::endl;
    return 0;
}
