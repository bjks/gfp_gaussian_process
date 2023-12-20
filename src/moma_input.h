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
#include "Gaussians.h"


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
    Eigen::VectorXi segment;

    // settings
    std::string noise_model;
    std::string cell_division_model;
    double fp_auto;

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

    // correlation function
    std::vector<Eigen::MatrixXd> correlation;
    Gaussian joint;
    

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

void build_cell_genealogy(std::vector<MOMAdata> &cell_vector){
    /*  
    * Assign respective pointers to parent, daughter1 and daughter2 for each cell
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
                    _file_log   << "(build_cell_genealogy) ERROR: Both daughter pointers are set, cell_id: " 
                                << cell_vector[j].cell_id << "\n";
                    _file_log << "-> daughter1 " << cell_vector[j].daughter1->cell_id << "\n";
                    _file_log << "-> daughter2 " << cell_vector[j].daughter2->cell_id << "\n";
                    throw std::invalid_argument("Invalid argument");
                }
            }
        }
    }
}

void print_cells(std::vector<MOMAdata> const &cell_vector){
    /* prints all cells */
    for (MOMAdata cell: cell_vector){
        _file_log << cell;
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
    std::string stri;
    for (size_t i = 0; i < str_vec.size(); ++i){
        stri = trim(str_vec[i], ' ');
        stri = trim(stri, '\t');
        stri = trim(stri, '\n');
        stri = trim(stri, '\v');
        stri = trim(stri, '\f');
        stri = trim(stri, '\r');
        header_indices.insert(std::pair<std::string, int>(stri, i)); 
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

void append_vec(Eigen::VectorXi &v, int elem){
    /*  
    * push_back alternative for non std vector (with resize() and size()), 
    * probaly slow and should only be used to read the csv 
    * and create vector with data with unknow length
    */
    v.conservativeResize(v.size()+1);
    v[v.size()-1] = elem;
}

double last_element(Eigen::VectorXd &v){
    return v[v.size()-1];
}

std::vector<MOMAdata> read_data(std::string filename, CSVconfig &config, std::string noise_model, std::string cell_division_model){
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
        _file_log << "(read_data) ERROR: (time_col) is not an column in input file: " << config.time_col << "\n";
        throw std::invalid_argument("Invalid argument");
        return data;
    }
    if (!header_indices.count(config.length_col)){
        _file_log << "(read_data) ERROR: (length_col) is not an column in input file: " << config.length_col << "\n";
        throw std::invalid_argument("Invalid argument");
        return data;
    }
    if (!header_indices.count(config.fp_col)){
        _file_log << "(read_data) ERROR: (fp_col) is not an column in input file: " << config.fp_col << "\n";
        throw std::invalid_argument("Invalid argument");
        return data;
    }
    if (!config.segment_col.empty()){
        if (!header_indices.count(config.segment_col)){
            _file_log << "(read_data) ERROR: (segment_col) is not an column in input file: " << config.segment_col << "\n";
            throw std::invalid_argument("Invalid argument");
            return data;
        }
    }
    if (!config.filter_col.empty()){
        if (!header_indices.count(config.filter_col)){
            _file_log << "(read_data) ERROR: (filter_col) is not an column in input file: " << config.filter_col << "\n";
            throw std::invalid_argument("Invalid argument");
            return data;
        }
    }
    for(size_t i=0; i<config.cell_tags.size(); ++i){
        if (!header_indices.count(config.cell_tags[i])){
            _file_log << "(read_data) ERROR: at least one of (cell_tags) is not an column in input file: " << config.cell_tags[i] << "\n";
            throw std::invalid_argument("Invalid argument");   
            return data;
        }
    }
    for(size_t i=0; i<config.parent_tags.size(); ++i){
        if (!header_indices.count(config.parent_tags[i])){
            _file_log << "(read_data) ERROR: at least one of (parent_tags) is not an column in input file: " << config.parent_tags[i] << "\n";
            throw std::invalid_argument("Invalid argument");
            return data;
        }
    }
    
    // Iterate through each line and split the content using the delimeter then assign the 
    std::string last_cell = "";
    std::string curr_cell;

    int last_idx = -1;
    long line_count = 1; // 1, since header line is already read
    while (getline(file, line)) {
        ++line_count;
        try{
            line_parts = split_string_at(line, config.delm);
            if (config.filter_col.empty() || string2bool(line_parts[header_indices[config.filter_col]])){

                // compose the cell id of the cells using the cell_tags
                curr_cell = get_cell_id(line_parts, header_indices, config.cell_tags);

                if (last_cell != curr_cell){
                    last_idx++;
                    MOMAdata next_cell;
                    // add new MOMAdata instance to vector 
                    data.push_back(next_cell); 

                    data[last_idx].cell_id = curr_cell;
                    // compose the cell id of the parent using the parent_tags
                    data[last_idx].parent_id = get_cell_id(line_parts, header_indices, config.parent_tags);
                    data[last_idx].noise_model = noise_model;
                    data[last_idx].cell_division_model = cell_division_model;
                    data[last_idx].fp_auto = config.fp_auto;
                }

                append_vec(data[last_idx].time,  stod_reject_nan(line_parts[header_indices[config.time_col]])/config.rescale_time);

                if (config.length_islog)
                    append_vec(data[last_idx].log_length,  stod_reject_nan(line_parts[header_indices[config.length_col]]) );
                else
                    append_vec(data[last_idx].log_length,  log(stod_reject_nan(line_parts[header_indices[config.length_col]])) );

                append_vec(data[last_idx].fp,  stod_reject_nan(line_parts[header_indices[config.fp_col]]) );

                if (config.segment_col.empty()){
                    append_vec(data[last_idx].segment, 0); //in case there is only one segment, assign a dummy segment index
                }
                else{
                    append_vec(data[last_idx].segment,  std::stoi(line_parts[header_indices[config.segment_col]]) );
                }

                /* ============ */
                last_cell = curr_cell;
            }
        }
        catch(std::exception &e){
            _file_log << "(read_data) ERROR: Line no." << line_count \
                        << " ["  \
                        << curr_cell << " , "  \
                        << line_parts[header_indices[config.time_col]] << " , "  \
                        << line_parts[header_indices[config.length_col]] << " , "  \
                        << line_parts[header_indices[config.fp_col]] << " , "  \
                        << line_parts[header_indices[config.segment_col]]  \
                        << "]" \
                        << " cannnot be processed (" << e.what() <<")" << std::endl;
            throw;
        }
    }
    file.close();
    _file_log << last_idx + 1 << " cells and " << line_count << " data points found in file " << filename << std::endl; 
    return data;
}


long count_data_points(std::vector<MOMAdata> const &cells){
    long ndata_points = 0;
    for(size_t i=0; i<cells.size(); ++i){
        ndata_points += cells[i].time.size();
    }
    return ndata_points;
}

std::vector<int> get_segment_indices(std::vector<MOMAdata> cells){
    /* 
    * Checks if the segment indices are consecutive and start at 0, 
    * Also returns the indices in order of occurence in the data set which deterines in which order the segments are run, 
    * (although this does not matter)
    */
    std::vector<int> segs;
    for(size_t i=0; i<cells.size(); ++i){
        for(size_t t=0; t<cells[i].time.size(); ++t){
            if( !(std::find(segs.begin(), segs.end(), cells[i].segment[t]) != segs.end()) ){
                segs.push_back(cells[i].segment[t]);
            }
        }
    }
    if (*std::min_element(segs.begin(), segs.end()) != 0){
        _file_log << "(get_segment_indices) ERROR: The segment indices do not start at 0:";
        for (size_t i=0; i<segs.size(); ++i){
            _file_log << " " << segs[i];
        }
        _file_log << "\n";
        throw std::invalid_argument("Invalid argument");
    }

    if (segs.size()-1 != *std::max_element(segs.begin(), segs.end())){
        _file_log << "(get_segment_indices) ERROR: The segment indices are not consecutive:";
        for (size_t i=0; i<segs.size(); ++i){
            _file_log << " " << segs[i];
        }
        _file_log << "\n";
        throw std::invalid_argument("Invalid argument");
    }
    return segs;
}

int get_segment_file_number(std::vector<int> segment_indices, int i){
    /* check if multi-segment mode is run, just for file naming */
    if (segment_indices.size()>1)
        return i;
    else
        return -1;
}

std::vector<MOMAdata> get_segment(std::vector<MOMAdata> cells, int segment){
    /* 
    * Returns a vector of those cells that are in the requested segment, 
    * note that this functio does not respect any pointers to other cells and the genealogy of pointers needs to be build afterwards
    */
    std::vector<MOMAdata> cells_in_segment;
    for(size_t i=0; i<cells.size(); ++i){
        MOMAdata cell;
        cell.cell_id                = cells[i].cell_id;
        cell.parent_id              = cells[i].parent_id;
        cell.noise_model            = cells[i].noise_model;
        cell.cell_division_model    = cells[i].cell_division_model;
        cell.fp_auto                = cells[i].fp_auto;

        for(size_t t=0; t<cells[i].time.size(); ++t){
            if(cells[i].segment[t] == segment){
                append_vec(cell.time, cells[i].time[t]);
                append_vec(cell.log_length, cells[i].log_length[t]);
                append_vec(cell.fp, cells[i].fp[t]);
                append_vec(cell.segment, cells[i].segment[t]);

                /* in case the prediction part is already run, keep those */
                if (cells[i].mean_forward.size()){ //note all "prediction" vectors are the same size always
                    cell.mean_forward.push_back(cells[i].mean_forward[t]);
                    cell.cov_forward.push_back(cells[i].cov_forward[t]);

                    cell.mean_backward.push_back(cells[i].mean_backward[t]);
                    cell.cov_backward.push_back(cells[i].cov_backward[t]);

                    cell.mean_prediction.push_back(cells[i].mean_prediction[t]);
                    cell.cov_prediction.push_back(cells[i].cov_prediction[t]);
                }
            }
        }
        /* only keep the cell if it has any data points in the segments */
        if (cell.time.size()){
            cells_in_segment.push_back(cell);
        }
    }
    return cells_in_segment;
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

double vec_mean(std::vector<double> v){
    /* returns mean of std::vector */
    return std::accumulate(v.begin(), v.end(), 0.0) / v.size();

}
double vec_var(std::vector<double> v){
    /* returns variance of std::vector */
    double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
    return sq_sum / v.size() - pow(vec_mean(v), 2);
}


void init_cells_f(std::vector<MOMAdata> &cells){
    /* 
    * Inititalizes the mean vector and the covariance matrix of the ROOT cells estimated from 
    * the data using the FIRST time point for x and fp 
    */

    // Estimate initial x, g, and lambda
    std::vector<double> x0;
    std::vector<double> g0;

    for(size_t i=0; i<cells.size(); ++i){
        if (cells[i].time.size()>1){
            x0.push_back(cells[i].log_length(0));
            g0.push_back(cells[i].fp(0));
        }
    }
    double mean_x0 = vec_mean(x0);
    double mean_g0 = vec_mean(g0);

    double var_x0 =  vec_var(x0);
    double var_g0 =  vec_var(g0);

    std::vector<MOMAdata *> roots = get_roots(cells);
    for(size_t i=0; i<roots.size(); ++i){
        roots[i]->mean_init_forward << mean_x0, mean_g0, 0, 0;
        roots[i]->cov_init_forward = Eigen::MatrixXd::Zero(4, 4);
        roots[i]->cov_init_forward(0,0) = var_x0;
        roots[i]->cov_init_forward(1,1) = var_g0;
    }
}

void init_cells_r(std::vector<MOMAdata> &cells){
    /* 
    * Inititalizes the mean vector and the covariance matrix of the LEAF cells estimated from 
    * the data using the LAST time point for x and fp 
    */

    // Estimate initial x, g, and lambda
    std::vector<double> x0;
    std::vector<double> g0;

    for(size_t i=0; i<cells.size(); ++i){
        if (cells[i].time.size()>1){
            x0.push_back(cells[i].log_length(cells[i].log_length.size()-1));
            g0.push_back(cells[i].fp(cells[i].fp.size()-1));
        }
    }
    double mean_x0 = vec_mean(x0);
    double mean_g0 = vec_mean(g0);

    double var_x0 = vec_var(x0);
    double var_g0 = vec_var(g0);

    std::vector<MOMAdata *> leafs = get_leafs(cells);
    for(size_t i=0; i<leafs.size(); ++i){
        leafs[i]->mean_init_backward << mean_x0, mean_g0, 0, 0;
        leafs[i]->cov_init_backward = Eigen::MatrixXd::Zero(4, 4);
        leafs[i]->cov_init_backward(0,0) = var_x0;
        leafs[i]->cov_init_backward(1,1) = var_g0;
    }
}


void init_cells(std::vector<MOMAdata> &cells){
    init_cells_f(cells);
    init_cells_r(cells);
}
