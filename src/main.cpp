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

    Config config("config.txt");

    std::vector<MOMAdata> cells =  getData(argv[1], 
                                        config.time_col,
                                        config.length_col,
                                        config.fp_col,
                                        config.delm);

    build_cell_genealogy(cells);
    print_cell_genealogy(cells);

    for (MOMAdata cell: cells){
        if (cell.parent !=NULL ){
            std::cout   << cell.cell_id << "\n"
                        << " \t    - " << cell.parent->daughter1->cell_id  << std::endl;
            if (cell.parent->daughter2 != NULL) 
                std::cout   << " \t or - " << cell.parent->daughter2->cell_id << std::endl;
        }
    }

    // pvector(cells[0].time);
    // pvector(cells[0].length);
    // pvector(cells[0].fp);

    std::cout << "Done." << std::endl;
    return 0;
}
