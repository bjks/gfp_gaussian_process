#include <fstream>
#include <iostream>
#include <sstream>  

#ifndef PARAMETERS_H
#define PARAMETERS_H

class Parameter{
/* 
* Single parameter, that can be 
    free (fixed=false, bound=false) 
    bound (bound=true, fixed=false) 
    fixed (fixed=true)

All of which contain an initial (double) and a name (string).

Non-fixed parameters have a step.

Bound parameters have upper/lower (double) which are the respective bounds
*/
public:
    bool fixed;
    bool bound;
    bool free;

    double init;
    bool set;

    double final;
    bool minimized;

    double step;
    double lower;
    double upper;

    std::string name;
    Parameter(std::string name_, double lower_ = -HUGE_VAL, double upper_ = HUGE_VAL){
        name = name_;
        lower = lower_;
        upper = upper_;

        minimized = false;
        set = false;

        fixed = false;
        bound = false;
        free = false;
    }

    void set_paramter(std::vector<std::string> parts){
        std::vector<std::string> val_split;
        val_split = split_string_at(parts[1], ",");
        for (size_t i=0; i<val_split.size(); ++i){
            val_split[i] = trim(val_split[i]);
        }

        try{
            if (val_split.size() == 4){
                init =  stod_reject_nan(val_split[0]);
                step =  stod_reject_nan(val_split[1]);
                lower =  stod_reject_nan(val_split[2]);
                upper =  stod_reject_nan(val_split[3]);
                bound = true;
                set = true;

            } else if (val_split.size() == 1){
                init =  stod_reject_nan(val_split[0]);
                fixed = true;
                set = true;

            }  else if (val_split.size() == 2){
                init =  stod_reject_nan(val_split[0]);
                step =  stod_reject_nan(val_split[1]);
                free = true;
                set = true;

            } else{
                throw std::invalid_argument("Invalide number of arguments");
            }
        }
        catch(std::exception &e){
            _file_log << "(set_paramter) ERROR: Parameter settings of '" << parts[0] << "' cannot be processed (" << e.what() <<")" << std::endl;
            throw;
        }
    }
};


class Parameter_set{
/*  Notation in Athos thesis:
    ---------------------------
Growth rate fluctualtions params:
    mean_lambda;     = \bar \lambda
    gamma_lambda;    = \gamma_\lambda
    var_lambda;      = \sigma_\lambda^2

gfp fluctuation params
    mean_q;          = \bar q
    gamma_q;         = \gamma_q
    var_q;           = \sigma_q^2

    beta;            = \beta

variance guess for length and gfp
    var_x;           = \sigma_x^2
    var_g;           = \sigma_g^2

cell division:
    var_dx;          = \sigma_{dx}^2
    var_dg;          = \sigma_{dg}^2

*/
protected:
    Parameter mean_lambda = Parameter("mean_lambda", 0.);
    Parameter gamma_lambda = Parameter("gamma_lambda", 0.);
    Parameter var_lambda = Parameter("var_lambda", 0.);

    Parameter mean_q = Parameter("mean_q", 0.);
    Parameter gamma_q = Parameter("gamma_q", 0.);
    Parameter var_q = Parameter("var_q", 0.);
    
    Parameter beta = Parameter("beta", 0.);

    Parameter var_x = Parameter("var_x", 0.);
    Parameter var_g = Parameter("var_g", 0.);

    Parameter var_dx = Parameter("var_dx", 0.);
    Parameter var_dg = Parameter("var_dg", 0.);

public:
    std::vector<Parameter> all;

    Parameter_set(std::string filename) {
        std::ifstream fin(filename);
        std::string line;
        std::vector<std::string> parts;

        if(fin) {
            // Overwrite defaults if in config file
            while (getline(fin, line)) {
                if (line[0] != '#' && line.size()){
                    parts = split_string_at(line, "=");

                    // remove whitespaces from the ends
                    parts[0] = trim(parts[0]);

                    if (parts[0] == "mean_lambda"){
                        mean_lambda.set_paramter(parts);
                    } else if (parts[0] == "gamma_lambda"){
                        gamma_lambda.set_paramter(parts);
                    } else if (parts[0] == "var_lambda"){
                        var_lambda.set_paramter(parts);
                    } else if (parts[0] == "mean_q"){
                        mean_q.set_paramter(parts);
                    } else if (parts[0] == "gamma_q"){
                        gamma_q.set_paramter(parts);
                    } else if (parts[0] == "var_q"){
                        var_q.set_paramter(parts);
                    } else if (parts[0] == "beta"){
                        beta.set_paramter(parts);
                    } else if (parts[0] == "var_x"){
                        var_x.set_paramter(parts);
                    } else if (parts[0] == "var_g"){
                        var_g.set_paramter(parts);
                    } else if (parts[0] == "var_dx"){
                        var_dx.set_paramter(parts);
                    } else if (parts[0] == "var_dg"){
                        var_dg.set_paramter(parts);
                    }

                }
            }
        }
        // create vector containing all paramters in well-defined order
        all = {mean_lambda, gamma_lambda, var_lambda, mean_q, gamma_q, var_q, beta, var_x, var_g, var_dx, var_dg};
    }

    friend std::ostream& operator<<(std::ostream& os, const Parameter_set& params);

    bool check_if_complete();
    bool has_nonfixed();
    void set_final(std::vector<double> vals);
    std::vector<double> get_final();
    std::vector<double> get_init();
    std::vector<int> non_fixed();
    
    void const to_csv(std::ostream &file);
    void const to_csv(std::string outfile, std::ios_base::openmode mode = std::ios_base::out);
};


bool Parameter_set::check_if_complete(){
    for (size_t i=0; i<all.size(); ++i){
        if (!all[i].set){
            _file_log << "(check_if_complete) ERROR: Parameter " << all[i].name << " not found in parameter file\n";
            throw std::invalid_argument("Invalide argument");
            return false;
        }
    }
    return true;
}

bool Parameter_set::has_nonfixed(){ 
    for (size_t i=0; i<all.size(); ++i){
        if (!all[i].fixed){
            return true;
        }
    }
    return false;
}

void const Parameter_set::to_csv(std::ostream &file){
    file << "no,name,type,init,step,lower_bound,upper_bound,final\n";
    for (size_t i=0; i<all.size(); ++i){
        file <<  i << "," ;
        if (all[i].fixed){
            file  << all[i].name << ","
                << "fixed"<< ","
                << all[i].init << ", , , ,";
        } else if (all[i].bound){
            file <<  all[i].name << ","
                << "bound" << ","
                << all[i].init << ","
                << all[i].step << ","
                << all[i].lower << ","
                << all[i].upper << ",";
        } else {
            file <<  all[i].name << ","
                << "free" << ","
                << all[i].init << ","
                << all[i].step << ", , ,";
        }
        if (all[i].minimized){
            file << all[i].final;
        }
    file << "\n";    
    }
}

void const Parameter_set::to_csv(const std::string outfile, std::ios_base::openmode mode){
    std::ofstream file(outfile, mode);
    file << "no,name,type,init,step,lower_bound,upper_bound,final\n";
    for (size_t i=0; i<all.size(); ++i){
        file <<  i << "," ;
        if (all[i].fixed){
            file  << all[i].name << ","
                << "fixed"<< ","
                << all[i].init << ", , , ,";
        } else if (all[i].bound){
            file <<  all[i].name << ","
                << "bound" << ","
                << all[i].init << ","
                << all[i].step << ","
                << all[i].lower << ","
                << all[i].upper << ",";
        } else {
            file <<  all[i].name << ","
                << "free" << ","
                << all[i].init << ","
                << all[i].step << ", , ,";
        }
        if (all[i].minimized){
            file << all[i].final;
        }
    file << "\n";    
    }
    file.close();
}


void Parameter_set::set_final(std::vector<double> vals){
    for (size_t i=0; i<all.size(); ++i){
        all[i].final = vals[i];
        all[i].minimized = true;
    }
    
}

std::vector<double> Parameter_set::get_final(){
    std::vector<double> vals;
    for (size_t i=0; i<all.size(); ++i){
        if (all[i].minimized){
            vals.push_back(all[i].final);
        }
        else{
            vals.push_back(all[i].init);
        }
    }
    return vals;
}

std::vector<double> Parameter_set::get_init(){
    std::vector<double> vals;
    for (size_t i=0; i<all.size(); ++i){
        vals.push_back(all[i].init);
    }
    return vals;
}

std::vector<int> Parameter_set::non_fixed(){
    std::vector<int> idx;
    for (size_t i=0; i<all.size(); ++i){
        if (!all[i].fixed){
            idx.push_back(i);
        }
    }
    return idx;
}

std::ostream& operator<<(std::ostream& os, const Parameter_set& params){
    std::vector<int> column_widths {4, 15, 8, 20, 20, 15, 15};

    os  << pad_str("No", column_widths[0]) << pad_str("Name", column_widths[1])<< pad_str("Type", column_widths[2])
        << pad_str("Init", column_widths[3]) << pad_str("Step", column_widths[4]) << pad_str("Bounds", column_widths[5]) << "\n";

    for(size_t i=0; i<accumulate(column_widths.begin(), column_widths.end(), 0); ++i){
        os << "_";
    }
    os << "\n";

    for (size_t i=0; i<params.all.size(); ++i){
        os <<  pad_str(std::to_string(i) +":", column_widths[0]);
        if (params.all[i].fixed){
            os  << pad_str(params.all[i].name, column_widths[1]) 
                << pad_str("(fixed)", column_widths[2]) 
                << pad_str(params.all[i].init, column_widths[3])
                << pad_str("", column_widths[4]+column_widths[5]+column_widths[6]);
        } else if (params.all[i].bound){
            os <<  pad_str(params.all[i].name, column_widths[1]) 
                << pad_str("(bound)", column_widths[2]) 
                << pad_str(params.all[i].init, column_widths[3])
                << pad_str(params.all[i].step , column_widths[4]) 
                << pad_str(params.all[i].lower , column_widths[5])
                << pad_str(params.all[i].upper , column_widths[6]);
        } else {
            os <<  pad_str(params.all[i].name, column_widths[1]) 
                << pad_str("(free)", column_widths[2]) 
                << pad_str(params.all[i].init, column_widths[3])
                << pad_str(params.all[i].step, column_widths[4])
                << pad_str("", column_widths[5]+column_widths[6]) ;
        }
        if (params.all[i].minimized && !params.all[i].fixed){
            os << " -> "<< params.all[i].final;
        }
        os << "\n";
    }
    return os;
}

#endif 
