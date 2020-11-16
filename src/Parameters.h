#include <fstream>
#include <iostream>

class Parameter{
/* 
* Single parameter, that can be 
    free (fixed=false, bound=false) 
    bound (bound=true, fixed=false) 
    fixed (fixed=true)

All of which contain an initial value (double) and a name (string).

Non-fixed parameters have a step.

Bound parameters have upper/lower (double) which are the respective bounds
*/
public:
    bool fixed;
    bool bound;

    double init;
    bool set = false;

    double value;
    bool miminized = false;

    double step;
    double lower;
    double upper;

    std::string name;

    void set_paramter(std::vector<std::string> parts){
        name = parts[0];
        
        std::vector<std::string> val_split;
        val_split = split_string_at(parts[1], ",");
        for (int i=0; i<val_split.size(); ++i){
            val_split[i] = trim(val_split[i]);
        }

        if (val_split.size() == 4){
            init =  std::stod(val_split[0]);
            step =  std::stod(val_split[1]);
            lower =  std::stod(val_split[2]);
            upper =  std::stod(val_split[3]);
            fixed = false;
            bound = true;
            set = true;

        } else if (val_split.size() == 1){
            init =  std::stod(val_split[0]);
            fixed = true;
            set = true;

        }  else if (val_split.size() == 2){
            init =  std::stod(val_split[0]);
            step =  std::stod(val_split[1]);
            fixed = false;
            bound = false;
            set = true;

        } else{
            std::cout << "Warning number of vaues in parameter file, paramter not set " << name << std::endl;
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
    Parameter mean_lambda;
    Parameter gamma_lambda;
    Parameter var_lambda;

    Parameter mean_q;
    Parameter gamma_q;
    Parameter var_q;
    
    Parameter beta;

    Parameter var_x;
    Parameter var_g;

    Parameter var_dx;
    Parameter var_dg;

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
    void set_value(std::vector<double> vals);
    std::vector<double> get_values();
};


void Parameter_set::set_value(std::vector<double> vals){
    for (int i=0; i<all.size(); ++i){
        all[i].value =vals[i];
        all[i].miminized = true;
    }
    
}

std::vector<double>  Parameter_set::get_values(){
    std::vector<double> vals;
    for (int i=0; i<all.size(); ++i){
        vals.push_back(all[i].value);
    }
    return vals;
}

std::string pad_str(std::string s, const size_t num, const char paddingChar = ' '){
    if(num > s.size())
        s.insert(s.end(), num - s.size(), paddingChar);
    return s;
}

std::ostream& operator<<(std::ostream& os, const Parameter_set& params){
    for (int i=0; i<params.all.size(); ++i){
        if (params.all[i].set){
            os <<  pad_str(std::to_string(i), 2) << ": ";
            if (params.all[i].fixed){
                os <<  pad_str(params.all[i].name, 15) << " (fixed) = " << params.all[i].init;
            } else if (params.all[i].bound){
                os << pad_str(params.all[i].name, 15) << " (bound) = " 
                    << params.all[i].init << ", step: " << params.all[i].step 
                    << ", bounds: (" << params.all[i].lower << ", " << params.all[i].upper << ")";
            } else {
                os << pad_str(params.all[i].name, 15) << " (free)  = " 
                    << params.all[i].init << ", step: " << params.all[i].step;
            }
            if (params.all[i].miminized && !params.all[i].fixed){
                os << " -> "<< params.all[i].value;
            }
        }
        else{
            os <<  pad_str(std::to_string(i), 2) << ": WARNING parameter not found, please see Parameters.h for info" << params.all[i].name;
        }
        os << "\n";

    }
    return os;
}

