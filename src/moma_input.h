#include <iostream>
#include <fstream>
#include <iterator>
#include <string>

#include <vector>
#include <map> 
#include <cmath>
#include <numeric> // for accumulate and inner_product

#include <Eigen/Core>
#include <Eigen/LU> 

#include "CSVconfig.h"

// ============================================================================= //
// MOMAdata CLASS
// ============================================================================= //

class MOMAdata{
    /*  
    * A class containing data from a MOMA-csv file (or similar) for a single cell 
    * that stores all calculated quanitites and its realtions to other cells in the data set
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

    int generation;

    // initial guess for mean and cov, to avoid recalculation
    Eigen::VectorXd mean_init_forward = Eigen::VectorXd::Zero(4);
    Eigen::MatrixXd cov_init_forward = Eigen::MatrixXd::Zero(4, 4);

    Eigen::VectorXd mean_init_backward = Eigen::VectorXd::Zero(4);
    Eigen::MatrixXd cov_init_backward = Eigen::MatrixXd::Zero(4, 4);

    // variables to for updating the current state
    Eigen::VectorXd mean = Eigen::VectorXd::Zero(4);
    Eigen::MatrixXd cov = Eigen::MatrixXd::Zero(4, 4);

    // predictions
    std::vector<Eigen::Vector4d> mean_forward;
    std::vector<Eigen::Matrix4d> cov_forward;

    std::vector<Eigen::Vector4d> mean_backward;
    std::vector<Eigen::Matrix4d> cov_backward;

    std::vector<Eigen::Vector4d> mean_prediction;
    std::vector<Eigen::Matrix4d> cov_prediction;

    // member functions
    bool is_leaf() const;
    bool is_root() const;

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
    os << cell.mean << "\n";
    os << cell.cov << "\n";

    return os;
}

bool MOMAdata :: is_leaf() const{
    /* returns true if cell is leaf in tree */
    return daughter1 == nullptr && daughter2 == nullptr;
}

bool MOMAdata :: is_root() const {
    /* returns true if cell is root in tree */
    return parent==nullptr;
}


// ============================================================================= //
// GENEALOGY
// ============================================================================= //

bool build_cell_genealogy(std::vector<MOMAdata> &cell_vector){
    /*  
    * Assign respective pointers to parent, daughter1 and daughter2 for each cell.
    * Returns true when each parent was assigned to max. 2 daughters. false otherwise
    */
    for(size_t k = 0; k < cell_vector.size(); ++k) {
        for(size_t j = 0; j < cell_vector.size(); ++j) {
            if( cell_vector[j].cell_id == cell_vector[k].parent_id ){
                //  Assign pointers to PARENT variable of the cell
                cell_vector[k].parent = &cell_vector[j];
                //  Assign pointers to CELL of the parent cell to 'free' pointer
                if (cell_vector[j].daughter1 == nullptr){
                    cell_vector[j].daughter1 = &cell_vector[k];
                }
                else if (cell_vector[j].daughter2 == nullptr){
                    cell_vector[j].daughter2 = &cell_vector[k];
                }
                else{
                    std::cerr   << "(build_cell_genealogy) ERROR: both daughter pointers are set!" 
                                << cell_vector[j].cell_id << "\n";
                    std::cerr << "-> daughter1 " << cell_vector[j].daughter1->cell_id << "\n";
                    std::cerr << "-> daughter2 " << cell_vector[j].daughter2->cell_id << "\n";
                    return false;   
                }
            }
        }
    }
    return true;
}

void print_cells(std::vector<MOMAdata> const &cell_vector){
    /* prints all cells */
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
    for(size_t i=0; i < cells.size(); ++i){
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
    for(size_t i=0; i < cells.size(); ++i){
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

// DOWN ------------------------------------------------------------------------ //
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

	       1            |
	     /   \          |
	    2     5         |
	  /   \     \       |
	 3     4     6      V

    * _________________________________________________ 
    */
    apply_down_tree_recr(params_vec, &cell, func);
}


// UP ------------------------------------------------------------------------ //
void apply_up_tree_recr(const std::vector<double> &params_vec, 
                        MOMAdata *cell, 
                        void (*func)(const std::vector<double> &, MOMAdata &))
                        {
    /*  
    * Recursive implementation that applies the function func to every cell in the genealogy
    * not meant to be called directly, see wrapper below
    */
    if (cell == nullptr)
        return;

    apply_down_tree_recr(params_vec, cell->daughter1, func);
    apply_down_tree_recr(params_vec, cell->daughter2, func);

    func(params_vec, *cell);

}

void apply_up_tree(const std::vector<double> &params_vec, 
                    MOMAdata &cell, 
                    void (*func)(const std::vector<double> &, MOMAdata &))
                    {
    /* applies the function func to the cell and the other cells in the genealogy
    * such that the daughter cells has already been accessed when the function is applied 
    * to the cell.
    * 
    * Example (number implies the order in which cell is accessed)
    * _________________________________________________ 

	       6            ^
	     /   \          |
	    3     5         |
	  /   \     \       |
	 1     2     4      |

    * _________________________________________________ 
    */
    apply_up_tree_recr(params_vec, &cell, func);
}



// ============================================================================= //
// READING CSV
// ============================================================================= //
std::string remove_last_decimal(std::string str){
    /* removes endings .0 .00 .000... of purely numeric strings */

    // check if only numeric chars in str
    for (size_t i = 0; i < str.size(); ++i){
        if (!isdigit(str[i]) && str[i] != '.')
            return str; 
    }

    // check if all characters after last '.' are 0s
    std::vector parts = split_string_at(str, ".");
    std::string last_part = parts[parts.size()-1];
    for (size_t i = 0; i < last_part.size(); ++i){
        if (last_part[i] != '0'){
            return str;
        }
    }
    return std::to_string(std::stoi(str));
}


std::string get_cell_id(std::vector<std::string> &str_vec, 
                            std::map<std::string, int> &header_indices,
                            std::vector<std::string> tags){
    /*  
    * Compose id of the cell by adding all elements in tags seperated by "."
    */
    std::string id =""; 
    for (size_t i=0; i<tags.size(); ++i){
        if (i>0)
            id += ".";
        id += remove_last_decimal(str_vec[header_indices[tags[i]]]);
        // need to get rid of decimal points, hence the double type cast
    }
    return id;
}


std::map<std::string, int> get_header_indices(std::vector<std::string> &str_vec){
    /*  
    * Create a map containing the header tags and the corresponding index
    */
    std::map<std::string, int> header_indices; 
    for (size_t i = 0; i < str_vec.size(); ++i){
        header_indices.insert(std::pair<std::string, int>(str_vec[i], i)); 
    }
    return header_indices;
}


void append_vec(Eigen::VectorXd &v, double elem){
    /*  
    * push_back alternative for non std vector (with resize() and size()), 
    * probaly slow and should only be used to read the csv 
    * and create vector with data with unknow length
    */
    v.conservativeResize(v.size()+1);
    v[v.size()-1] = elem;
}


std::vector<MOMAdata> get_data(std::string filename, CSVconfig &config){
    /* 
    * Parses csv file line by line and returns the data as a vector of MOMAdata instances.
    * Returns data as vector of MOMAdata instances. Pointers for genealogy are not set yet!
    */
    std::ifstream file(filename);
    
    std::vector<std::string> line_parts;
    std::string line;
    std::vector<MOMAdata> data;

    // read the header and assign an index to every entry, such that we can 'index' with a string
    getline(file, line);
    line_parts = split_string_at(line, config.delm);
    std::map<std::string, int> header_indices = get_header_indices(line_parts);

    // check if the columns that are set actually exist in header 
    if (!header_indices.count(config.time_col)){
        std::cerr << "(get_data) ERROR: (time_col) is not an column in input file: " << config.time_col << "\n";
        return data;
    }
    if (!header_indices.count(config.length_col)){
        std::cerr << "(get_data) ERROR: (length_col) is not an column in input file: " << config.length_col << "\n";
        return data;
    }
    if (!header_indices.count(config.fp_col)){
        std::cerr << "(get_data) ERROR: (fp_col) is not an column in input file: " << config.fp_col << "\n";
        return data;
    }
    for(size_t i=0; i<config.cell_tags.size(); ++i){
        if (!header_indices.count(config.cell_tags[i])){
            std::cerr << "(get_data) ERROR: at least one of (cell_tags) is not an column in input file: " << config.cell_tags[i] << "\n";
            return data;
        }
    }
    for(size_t i=0; i<config.parent_tags.size(); ++i){
        if (!header_indices.count(config.parent_tags[i])){
            std::cerr << "(get_data) ERROR: at least one of (parent_tags) is not an column in input file: " << config.parent_tags[i] << "\n";
            return data;
        }
    }
    
    // Iterate through each line and split the content using the delimeter then assign the 
    std::string last_cell = "";
    std::string curr_cell;

    int last_idx = -1;
    long line_count = 0;
    long data_point_cell;
    while (getline(file, line)) {
        ++line_count;
        line_parts = split_string_at(line, config.delm);
        // take lines only if end_type==div or header_indices "end_type" is not in header_indices
        if (header_indices.count("end_type") == 0 || line_parts[header_indices["end_type"]] == "div" ){
            // compose the cell id of the cells using the cell_tags
            curr_cell = get_cell_id(line_parts, header_indices, config.cell_tags);

            if (last_cell != curr_cell){
                data_point_cell = 0;
                last_idx++;
                MOMAdata next_cell;
                // add new MOMAdata instance to vector 
                data.push_back(next_cell); 

                data[last_idx].cell_id = curr_cell;
                // compose the cell id of the parent using the parent_tags
                data[last_idx].parent_id = get_cell_id(line_parts, header_indices, config.parent_tags);
            }

            append_vec(data[last_idx].time,  std::stod(line_parts[header_indices[config.time_col]])/config.rescale_time);

            if (config.length_islog)
                append_vec(data[last_idx].log_length,  std::stod(line_parts[header_indices[config.length_col]]) );
            else
                append_vec(data[last_idx].log_length,  log(std::stod(line_parts[header_indices[config.length_col]])) );

            append_vec(data[last_idx].fp,  std::stod(line_parts[header_indices[config.fp_col]]) );
            last_cell = curr_cell;
        }
    }
    file.close();
    std::cout << last_idx + 1 << " cells and " << line_count << " data points found in file " << filename << std::endl; 
    return data;
}


long count_data_points(std::vector<MOMAdata> const &cells){
    long ndata_points = 0;
    for(size_t i=0; i<cells.size(); ++i){
        ndata_points += cells[i].time.size();
    }
    return ndata_points;
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

double lin_fit_intercept(Eigen::VectorXd x, Eigen::VectorXd y) {
    /* returns intercept of linear regression */
    double s_x  = x.sum();
    double s_y  = y.sum();
    return (x.dot(x)*s_y - s_x*x.dot(y)) / (x.dot(x) * x.size() - s_x *s_x);
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

double estimate_lambda(MOMAdata &cell){
    /* Straigt-forward estimation of the growth rate taking the last and the first data point of cell */
    size_t t_last = cell.time.size()-1;
    double l = (cell.log_length(t_last) - cell.log_length(0)) / (cell.time(t_last) - cell.time(0));
    return l;
}

double estimate_q(MOMAdata &cell, double lambda_est){
    /* Straigt-forward estimation of q taking the last and the first data point of cell, while assuming beta=0 */
    size_t t_last = cell.time.size()-1;
    double dg = cell.fp(t_last) - cell.fp(0);
    double dv = exp(cell.log_length(t_last)) - exp(cell.log_length(0));
    double q = dg*lambda_est/dv;
    return q;
}

double estimate_q_stationary(MOMAdata &cell){
    /* For non growing cells the estimation of q via lambda is critical, hence this version for lambda -> 0 */
    size_t t_last = cell.time.size()-1;
    double dg = cell.fp(t_last) - cell.fp(0);
    double v0 = exp(cell.log_length(0));
    double q = dg/(v0* (cell.time(t_last) - cell.time(0)) ) ;
    return q;
}

void init_cells(std::vector<MOMAdata> &cells, bool stationary){
    /* 
    * Inititalizes the mean vector and the covariance matrix of the ROOT cells estimated from 
    * the data using the FIRST time point for x and fp 
    */

    // Estimate initial x, g, and lambda
    std::vector<double> x0;
    std::vector<double> g0;
    std::vector<double> l0;

    for(size_t i=0; i<cells.size(); ++i){
        if (cells[i].time.size()>1){
            x0.push_back(cells[i].log_length(0));
            g0.push_back(cells[i].fp(0));
            l0.push_back(estimate_lambda(cells[i]));
        }
    }
    double mean_x0 = vec_mean(x0);
    double mean_g0 = vec_mean(g0);
    double mean_l0 = vec_mean(l0);
    
    // Estimate initial q, which needs some guess for lambda
    std::vector<double> q0;
    for(size_t i=0; i<cells.size(); ++i){
        if (cells[i].time.size()>1){
            if (stationary){
                q0.push_back(estimate_q_stationary(cells[i]));
            }
            else{
                q0.push_back(estimate_q(cells[i], mean_l0));
            }
        }
    }
    double mean_q0 = vec_mean(q0);

    double var_x0 =  vec_var(x0);
    double var_g0 =  vec_var(g0);
    double var_l0 =  vec_var(l0);
    double var_q0 =  vec_var(q0);

    std::vector<MOMAdata *> roots = get_roots(cells);
    for(size_t i=0; i<roots.size(); ++i){
        roots[i]->mean_init_forward << mean_x0, mean_g0, mean_l0, mean_q0;

        roots[i]->cov_init_forward = Eigen::MatrixXd::Zero(4, 4);
        roots[i]->cov_init_forward(0,0) = var_x0;
        roots[i]->cov_init_forward(1,1) = var_g0;
        roots[i]->cov_init_forward(2,2) = var_l0;
        roots[i]->cov_init_forward(3,3) = var_q0;
    }
}


void init_cells_r(std::vector<MOMAdata> &cells, bool stationary){
    /* 
    * Inititalizes the mean vector and the covariance matrix of the LEAF cells estimated from 
    * the data using the LAST time point for x and fp 
    */

    // Estimate initial x, g, and lambda
    std::vector<double> x0;
    std::vector<double> g0;
    std::vector<double> l0;

    for(size_t i=0; i<cells.size(); ++i){
        if (cells[i].time.size()>1){
            x0.push_back(cells[i].log_length(cells[i].log_length.size()-1));
            g0.push_back(cells[i].fp(cells[i].fp.size()-1));
            l0.push_back(estimate_lambda(cells[i]));
        }
    }
    double mean_x0 = vec_mean(x0);
    double mean_g0 = vec_mean(g0);
    double mean_l0 = vec_mean(l0);
    
    // Estimate initial q, which needs some guess for lambda
    std::vector<double> q0;
    for(size_t i=0; i<cells.size(); ++i){
        if (cells[i].time.size()>1){
            if (stationary){
                q0.push_back(estimate_q_stationary(cells[i]));
            }
            else{
                q0.push_back(estimate_q(cells[i], mean_l0));
            }
        }
    }
    double mean_q0 = vec_mean(q0);

    double var_x0 = vec_var(x0);
    double var_g0 = vec_var(g0);
    double var_l0 = vec_var(l0);
    double var_q0 = vec_var(q0);

    std::vector<MOMAdata *> leafs = get_leafs(cells);
    for(size_t i=0; i<leafs.size(); ++i){
        leafs[i]->mean_init_backward << mean_x0, mean_g0, -mean_l0, -mean_q0;

        leafs[i]->cov_init_backward = Eigen::MatrixXd::Zero(4, 4);

        leafs[i]->cov_init_backward(0,0) = var_x0;
        leafs[i]->cov_init_backward(1,1) = var_g0;
        leafs[i]->cov_init_backward(2,2) = var_l0;
        leafs[i]->cov_init_backward(3,3) = var_q0;
    }
}


void init_cells(std::vector<MOMAdata> &cells, Eigen::VectorXd mean, Eigen::MatrixXd cov){
    /* 
    * Inititalizes the mean vector and the covariance matrix of the root cells with
    * pre-defined values
    */
    for(size_t i=0; i<cells.size(); ++i){
        cells[i].mean_init_forward = Eigen::VectorXd::Zero(4);
        cells[i].cov_init_forward = Eigen::MatrixXd::Zero(4, 4);
    }

    std::vector<MOMAdata *> roots = get_roots(cells);
    for(size_t i=0; i<roots.size(); ++i){
        roots[i]->mean_init_forward = mean;
        roots[i]->cov_init_forward = cov;
    }
}
