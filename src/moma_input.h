#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <string>
#include <algorithm>
#include <boost/algorithm/string.hpp>


class MOMAdata{
    /*  
    * A class containing data from a MOMA-csv file
    */
public:
    // IDs (eg '20150624.0.1.5') of related cells saved as strings
    std::string cell_id;
    std::string parent_id;

    // Pointer to other instances of the class representing the genealogy
    MOMAdata *parent;
    MOMAdata *daughter1;
    MOMAdata *daughter2;

    // Time dependent quantities (and time) of the cell
    std::vector<double> time;
    std::vector<double> length;
    std::vector<double> fp;

};


std::string get_parent_id(std::vector<std::string> &str_vec, std::map<std::string, int> &header_indices){
    /*  
    * Compose parent_id of the cell
    */
    std::string parent_id = str_vec[header_indices["date"]] + "." + 
                            str_vec[header_indices["pos"]]+ "." + 
                            str_vec[header_indices["gl"]] + "." + 
                            std::to_string(std::stoi(str_vec[header_indices["parent_id"]])); 
                            // need to get rid of decimals in "parent_id", hence the double type cast
    return parent_id;
}


std::map<std::string, int> get_header_indices(std::vector<std::string> &str_vec){
    /*  
    * Create a map containing the header tags and the corresponding index
    */
    std::map<std::string, int> header_indices; 
    for (int i = 0; i < str_vec.size(); ++i){
        header_indices.insert(std::pair<std::string, int>(str_vec[i], i)); 
    }
    return header_indices;
}


void build_cell_genealogy(std::vector<MOMAdata> &cell_vector){
    /*  
    * Assign respective pointers to parent, daughter1 and daughter2 for each cell
    */
    for(long k = 0; k < cell_vector.size(); ++k) {
        for(long j = 0; j < cell_vector.size(); ++j) {

            if( cell_vector[j].cell_id == cell_vector[k].parent_id ){
                //  Assign pointers to PARENT variable of the cell
                cell_vector[k].parent = &cell_vector[j];
                //  Assign pointers to CELL of the parent cell to 'free' pointer
                if (cell_vector[j].daughter1 == NULL)
                    cell_vector[j].daughter1 = &cell_vector[k];
                else if (cell_vector[j].daughter2 == NULL)
                    cell_vector[j].daughter2 = &cell_vector[k];
            }

        }
    }
}


void print_cell_genealogy(std::vector<MOMAdata> &cell_vector){
    for (MOMAdata cell: cell_vector){
        if (cell.parent !=NULL)
            std::cout << cell.cell_id << " \t <- parent: " << cell.parent->cell_id << std::endl;
        else
            std::cout << cell.cell_id << std::endl;

        if (cell.daughter1 !=NULL)
            std::cout << "\t -> daughter 1: " << cell.daughter1->cell_id << std::endl;     
        if (cell.daughter2 !=NULL)
            std::cout << "\t -> daughter 2: " << cell.daughter2->cell_id << std::endl;     
    }
}


std::vector<MOMAdata> getData(std::string filename,
                            std::string time_col, 
                            std::string length_col, 
                            std::string fp_col, 
                            std::string delm){
    /*  
    * Parses through csv file line by line and returns the data as a vector of MOMAdata instances
    */
    std::ifstream file(filename);
    
    std::vector<std::string> vec;
    std::string line;

    // read the header and assign an index to every entry, such that we can 'index' with a string
    getline(file, line);
    boost::algorithm::split(vec, line, boost::is_any_of(delm));
    std::map<std::string, int> header_indices = get_header_indices(vec);

    // Iterate through each line and split the content using the delimeter then assign the 
    std::string last_cell = "";
    std::string curr_cell;

    std::vector<MOMAdata> data;
    int last_idx = -1;
    while (getline(file, line)) {
        boost::algorithm::split(vec, line, boost::is_any_of(delm));
        curr_cell = vec[header_indices["cell"]];

        if (last_cell != curr_cell){
            last_idx++;
            MOMAdata next_cell;
            // add new MOMAdata instance to vector 
            data.push_back(next_cell); 

            data[last_idx].cell_id = curr_cell;
            data[last_idx].parent_id = get_parent_id(vec, header_indices);
        }
        data[last_idx].time.push_back(std::stod(vec[header_indices[time_col]]));
        data[last_idx].length.push_back(std::stod(vec[header_indices[length_col]]));
        data[last_idx].fp.push_back(std::stod(vec[header_indices[fp_col]]));

        last_cell = curr_cell;
    }
    file.close();
    return data;
}

