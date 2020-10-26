#include "read_csv.h"
#include "moma_input.h"

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

    // std::vector<MOMAdata> cells =  data.getData("myframes_asc662_lactose.csv");

    std::cout << "Input file " << argv[1] << std::endl;
    std::vector<MOMAdata> cells = getData(argv[1]);
    build_cell_genealogy(cells);
    print_cell_genealogy(cells);

    for (MOMAdata cell: cells){
        if (cell.parent !=NULL ){
            std::cout   << cell.cell_id << "\n"
                        << " \t    - " << cell.parent->daughter1->cell_id  << std::endl;
            if (cell.parent->daughter2 != NULL) 
                std::cout   << " \t or - " << cell.parent->daughter2->cell_id 
                            << std::endl;
        }
    }
    std::cout << "Done." << std::endl;
    return 0;
}
