#include "moma_input.h"
#include "read_config.h"

#include <iostream> 
#include <iterator> 
#include <map> 

void pvector(std::vector <std::string> const &a) {
    for(int i=0; i < a.size(); i++){
        std::cout << a[i] << ' ';
    }
    std::cout << std::endl;
}

void pvector(std::vector <double> const &a) {
    for(int i=0; i < a.size(); i++){
        std::cout << a[i] << ' ';
    }
    std::cout << std::endl;
}


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

    // get the "tree" starting from cell[0]
    std::vector<std::vector<MOMAdata *> > cell_paths = get_genealogy(&cells[0]);

    std::cout << std::endl;
    for (std::vector<MOMAdata *> path : cell_paths){
        for (MOMAdata * entry : path){
            std::cout << " -> " << entry->cell_id ;
        }
        std::cout << std::endl;
    }


    // print the data of the first cell
    pvector(cells[0].time);
    pvector(cells[0].length);
    pvector(cells[0].fp);


    std::cout << "Done." << std::endl;
    return 0;
}
