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

    // Read data
    std::string infile = argv[1];
    std::vector<MOMAdata> cells =  getData(infile, 
                                        config.time_col,
                                        config.length_col,
                                        config.fp_col,
                                        config.delm);
    build_cell_genealogy(cells);

    // some output to see cell relations
    print_cell_genealogy(cells);

    for (MOMAdata cell: cells){
        if (cell.parent !=NULL ){
            std::cout   << cell.cell_id << "\n"
                        << " \t    - " << cell.parent->daughter1->cell_id  << std::endl;
            if (cell.parent->daughter2 != NULL) 
                std::cout   << " \t or - " << cell.parent->daughter2->cell_id << std::endl;
        }
    }
    /*
    20150624.0.1.2
 	    - 20150624.0.1.2
 	 or - 20150624.0.1.4
    20150624.0.1.3
            - 20150624.0.1.3
        or - 20150624.0.1.5
    20150624.0.1.4
            - 20150624.0.1.2
        or - 20150624.0.1.4
    20150624.0.1.5
            - 20150624.0.1.3
        or - 20150624.0.1.5
    20150624.0.1.6
            - 20150624.0.1.6
    */

    // print the data of the first cell
    pvector(cells[0].time);
    pvector(cells[0].length);
    pvector(cells[0].fp);
    /*  
    7200 
    1.7615 
    4146.41 
    */

    std::cout << "Done." << std::endl;
    return 0;
}
