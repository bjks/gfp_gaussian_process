/* genealogy testing */

// example function to illustrate how genealogy realtions are handled here
void set_generation(const std::vector<double> &params_vec, MOMAdata &cell){
    if (cell.parent != nullptr){
        cell.generation = cell.parent->generation + 1;
    } else{
        cell.generation = 0;
    }
}

void set_generation_r(const std::vector<double> &params_vec, MOMAdata &cell){
    if (cell.daughter1 != nullptr){
        cell.generation = cell.daughter1->generation + 1;
    } else{
        cell.generation = 0;
    }
}


void print_generation_tree(const std::vector<double> &params_vec, MOMAdata &cell, std::string direction = "down"){

    // apply the set_generation function to root cell followed by the first generation etc...
    if (direction == "up" )
        apply_up_tree(params_vec, cell, set_generation_r);   
    else
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
