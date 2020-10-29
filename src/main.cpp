#include "read_config.h"
#include "moma_input.h"
#include "utils.h"
#include "likelihood.h"



#include <iostream> 
#include <iterator> 
#include <map> 


int main(int argc, char** argv){

    std::string infile = argv[1];
    
    if(! std::__fs::filesystem::exists(infile)){
        std::cout << "File " << infile << " not found! \nQuit" << std::endl;
        return 0;
    }

    Config config("config.txt");

    // Read data
    std::vector<MOMAdata> cells =  getData(infile, 
                                        config.time_col,
                                        config.length_col,
                                        config.fp_col,
                                        config.delm);
    build_cell_genealogy(cells);
    print_related_cells(cells);

    std::cout << cells[0].is_root() << std::endl;


    apply_down_tree(cells[0], set_generation);

    // get the "tree" starting from cell[0]
    std::vector<std::vector<MOMAdata *> > cell_paths = get_genealogy_paths(cells[0]);

    std::cout << std::endl;
    for (std::vector<MOMAdata *> path : cell_paths){
        for (MOMAdata * cell : path){
            std::cout << " -> " << cell->cell_id << " generation: " << cell->generation  ;
        }
        std::cout << std::endl;
    }
    
    /*
    -> 20150624.0.1.0 generation: 0 -> 20150624.0.1.2 generation: 1 -> 20150624.0.1.6 generation: 2
    -> 20150624.0.1.0 generation: 0 -> 20150624.0.1.4 generation: 1
    */


    // print the data of the first cell
    pvector(cells[0].time);
    pvector(cells[0].length);
    pvector(cells[0].fp);



    std::cout << "Done." << std::endl;
    return 0;
}
