#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <string>
#include <map> 
#include <math.h> 
#include <cmath>
#include <numeric>

#include "Parameters.h"

#include <Eigen/Core>
#include <Eigen/LU> 

// ============================================================================= //
// MOMAdata CLASS
// ============================================================================= //

class MOMAdata{
    /*  
    * A class containing data from a MOMA-csv file for a single cell
    */
public:
    // IDs (eg '20150624.0.1.5') of related cells saved as strings
    std::string cell_id;
    std::string parent_id;

    // Pointer to other instances of the class representing the genealogy
    MOMAdata *parent = nullptr;
    MOMAdata *daughter1 = nullptr;
    MOMAdata *daughter2 = nullptr;

    // Time dependent quantities (and time) of the cell
    // stores in eigen vectors to enable lin algebra functions
    Eigen::VectorXd time;
    Eigen::VectorXd log_length;
    Eigen::VectorXd fp;

    int generation; // to be deleted later 

    // member functions
    bool is_leaf() const;
    bool is_root() const;

    // variables to be calculated
    double likelihood;
    Eigen::VectorXd mean = Eigen::VectorXd::Zero(4);
    Eigen::MatrixXd cov = Eigen::MatrixXd::Zero(4, 4);


    friend std::ostream& operator<<(std::ostream& os, const MOMAdata& cell);
};


std::ostream& operator<<(std::ostream& os, const MOMAdata& cell){
    /*
    example output:
    ---------------
    20150624.0.1.0
        -> daughter 1: 20150624.0.1.2
        -> daughter 2: 20150624.0.1.4
    20150624.0.1.2 	 <- parent: 20150624.0.1.0
        -> daughter 1: 20150624.0.1.6
    20150624.0.1.3 	 <- parent: 20150624.0.1.1
    */
    
    os << cell.cell_id;
    if (cell.parent != nullptr)
        os << " \t <- parent: " << cell.parent->cell_id;
    os << "\n";

    if (cell.daughter1 !=nullptr)
        os << "\t \\_ daughter 1: " << cell.daughter1->cell_id << "\n";     
    if (cell.daughter2 !=nullptr)
        os << "\t \\_ daughter 2: " << cell.daughter2->cell_id << "\n";   
    return os;
}

bool MOMAdata :: is_leaf() const{
    return daughter1 == nullptr && daughter2 == nullptr;
}

bool MOMAdata :: is_root() const {
    return parent==nullptr;
}


// ============================================================================= //
// GENEALOGY
// ============================================================================= //

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
                if (cell_vector[j].daughter1 == nullptr)
                    cell_vector[j].daughter1 = &cell_vector[k];
                else if (cell_vector[j].daughter2 == nullptr)
                    cell_vector[j].daughter2 = &cell_vector[k];
                else
                    std::cout << "(build_cell_genealogy) Warning: both daughter pointers are set!" << std::endl;
            }
        }
    }
}

void print_cells(std::vector<MOMAdata> const &cell_vector){
    for (MOMAdata cell: cell_vector){
        std::cout << cell;
    }
}

// ----------------------------------------------------------------------------- //
// genealogy operations
// ----------------------------------------------------------------------------- //
std::vector<MOMAdata *> get_leafs(std::vector<MOMAdata > &cells){
    /*
    * returns vector pointers to MOMAdata cells 
    * each pointer points to a leaf of the cell tree
    */
    std::vector<MOMAdata *> leafs;
    for(int i=0; i < cells.size(); ++i){
        if (cells[i].is_leaf()){
            leafs.push_back(&cells[i]);
        }
    }
    return leafs;
}

std::vector<MOMAdata *> get_roots(std::vector<MOMAdata > &cells){
    /*
    * returns vector pointers to MOMAdata cells 
    * each pointer points to a root of the cell trees
    */
   
    std::vector<MOMAdata *> roots;
    for(int i=0; i < cells.size(); ++i){
        if (cells[i].is_root()){
            roots.push_back(&cells[i]);
        }
    }
    return roots;
}



// ----------------------------------------------------------------------------- //
// recursive path finding 
// ----------------------------------------------------------------------------- //

void get_genealogy_paths_recr(MOMAdata *cell, 
                                std::vector<MOMAdata *> &current_path, 
                                std::vector<std::vector<MOMAdata *> > &paths){
    /*
    * recursive function called by get_genealogy_paths which wrappes this one
    * not meant to be called directly, see wrapper below
    */
    if (cell == nullptr)
        return;

    current_path.push_back(cell);

    if (cell->is_leaf()){
        paths.push_back(current_path);
    } else{  
        get_genealogy_paths_recr(cell->daughter1, current_path, paths);
        get_genealogy_paths_recr(cell->daughter2, current_path, paths);
    }

    current_path.pop_back();
}

std::vector<std::vector<MOMAdata *> > get_genealogy_paths(MOMAdata &cell){
    /*
    * returns vector of vectors of pointers to MOMAdata cells 
    * each vector is a path from the given cell to one of its leafs
    */
    std::vector<MOMAdata *> current_path;
    std::vector<std::vector<MOMAdata *> > paths;

    get_genealogy_paths_recr(&cell, current_path, paths);
    return paths;
}

// ----------------------------------------------------------------------------- //
// recursive "looping"
// ----------------------------------------------------------------------------- //

void apply_down_tree_recr(const std::vector<double> &params_vec, 
                        MOMAdata *cell, 
                        void (*func)(const std::vector<double> &, MOMAdata &))
                        {
    /*  
    * Recursive implementation that applies the function func to every cell in the genealogy
    * not meant to be called directly, see wrapper below
    */
    if (cell == nullptr)
        return;
    func(params_vec, *cell);

    apply_down_tree_recr(params_vec, cell->daughter1, func);
    apply_down_tree_recr(params_vec, cell->daughter2, func);
}

void apply_down_tree(const std::vector<double> &params_vec, 
                    MOMAdata &cell, 
                    void (*func)(const std::vector<double> &, MOMAdata &))
                    {
    /* applies the function func to the cell cell and the other cells in the genealogy
    * such that the parent cell has already been accessed when the function is applied 
    * to the cell.
    * 
    * Example (number implies the order in which)
    * _________________________________________________ 

	       1
	     /   \
	    2     5
	  /   \     \
	 3     4     6

    * _________________________________________________ 
    */
    apply_down_tree_recr(params_vec, &cell, func);
}



// ============================================================================= //
// READING CSV
// ============================================================================= //
std::string get_parent_id(std::vector<std::string> &str_vec, 
                            std::map<std::string, int> &header_indices){
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


void append_vec(Eigen::VectorXd &v, double elem){
    /*  
    * push_back alternative for non std vector (with resize() and size()), probaly slow and should only be used to read the csv 
    * and create vector with data with unknow length
    */
    v.conservativeResize(v.size()+1);
    v[v.size()-1] = elem;
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
    
    std::vector<std::string> line_parts;
    std::string line;

    // read the header and assign an index to every entry, such that we can 'index' with a string
    getline(file, line);
    line_parts = split_string_at(line, delm);
    std::map<std::string, int> header_indices = get_header_indices(line_parts);

    // Iterate through each line and split the content using the delimeter then assign the 
    std::string last_cell = "";
    std::string curr_cell;

    std::vector<MOMAdata> data;
    int last_idx = -1;
    while (getline(file, line)) {
        line_parts = split_string_at(line, delm);
        if (line_parts[header_indices["end_type"]] == "div"){
            curr_cell = line_parts[header_indices["cell"]];

            if (last_cell != curr_cell){
                last_idx++;
                MOMAdata next_cell;
                // add new MOMAdata instance to vector 
                data.push_back(next_cell); 

                data[last_idx].cell_id = curr_cell;
                data[last_idx].parent_id = get_parent_id(line_parts, header_indices);
            }

            append_vec(data[last_idx].time,  std::stod(line_parts[header_indices[time_col]]) );
            append_vec(data[last_idx].log_length,  log(std::stod(line_parts[header_indices[length_col]]) ));
            append_vec(data[last_idx].fp,  std::stod(line_parts[header_indices[fp_col]]) );
            last_cell = curr_cell;
        }
    }
    file.close();
    std::cout << last_idx + 1 << " cells found in file " << filename << std::endl; 
    return data;
}

// ============================================================================= //
// MEAN/COV
// ============================================================================= //

Eigen::MatrixXd cov(Eigen::MatrixXd m){
    /*
    takes a matrix with rowise data, ie
        m = x1, x2, ...
            y1, y2, ...
            ...
    and calc covariance matrix between x,y,...
    */
    Eigen::MatrixXd cov;
    cov = m;
    Eigen::MatrixXd ones = Eigen::MatrixXd::Constant(1,m.cols(), 1);

    for(int i=0; i < m.rows(); ++i ){
        cov.row(i) -= ones * m.row(i).mean();
    } 
    return (cov * cov.transpose()) / (m.cols() - 1 );
}

Eigen::VectorXd row_mean(Eigen::MatrixXd m){
    /*
    takes a matrix with rowise data, ie
        m = x1, x2, ...
            y1, y2, ...
            ...
    and calc mean for x,y,...
    */
    Eigen::VectorXd mean(m.rows());
    for(int i=0; i < m.rows(); ++i ){
        mean(i) = m.row(i).mean();
    } 
    return mean;
}

// ============================================================================= //
// MEAN/COV INIT 
// ============================================================================= //

double lin_fit_slope(Eigen::VectorXd x, Eigen::VectorXd y) {
    /* returns slope of linear regression */
    double s_x  = x.sum();
    double s_y  = y.sum();
    return (x.size() * x.dot(y) - s_x * s_y) / (x.size() * x.dot(x) - s_x * s_x);
}

double vec_mean(std::vector<double> v){
    /* returns mean of std::vector */
    return std::accumulate(v.begin(), v.end(), 0.0) / v.size();

}
double vec_var(std::vector<double> v){
    /* returns variance of std::vector */
    double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
    return sq_sum / v.size() - pow(vec_mean(v), 2);
}

void init_cells(std::vector<MOMAdata> &cells){
    /* 
    * Inititalizes the mean vector and the covariance matrix of the root cells 
    */
    std::vector<double> x0;
    std::vector<double> g0;
    std::vector<double> l0;
    std::vector<double> q0;

    for(int i=0; i<cells.size(); ++i){
        if(cells[i].time.size()>2){
            x0.push_back(cells[i].log_length(0));
            g0.push_back(cells[i].fp(0));
            l0.push_back(lin_fit_slope(cells[i].time.head(3), cells[i].log_length.head(3)));
            q0.push_back(lin_fit_slope(cells[i].time.head(3), cells[i].fp.head(3)));
        }
    }

    std::vector<MOMAdata *> roots = get_roots(cells);
    for(int i=0; i<roots.size(); ++i){
        roots[i]->mean << vec_mean(x0),vec_mean(g0),vec_mean(l0),vec_mean(q0);
        roots[i]->cov(0,0) = vec_var(x0);
        roots[i]->cov(1,1) = vec_var(g0);
        roots[i]->cov(2,2) = vec_var(l0);
        roots[i]->cov(3,3) = vec_var(q0);
    }
}

