#include <fstream>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/numeric/ublas/vector.hpp>

class Parameter{
/* 
* Single parameter
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
        boost::algorithm::split(val_split, parts[1], boost::is_any_of(","));

        for (int i=0; i<val_split.size(); ++i){
            boost::algorithm::trim(val_split[i]);
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
            bound = true;
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
    mean_lambda;     = \bar \lambda
    gamma_lambda;    = \gamma_\lambda
    var_lambda;      = \sigma_\lambda^2

    mean_q;          = \bar q
    gamma_q;         = \gamma_q
    var_q;           = \sigma_q^2

    beta;            = \beta

    var_x;           = \sigma_x^2
    var_g;           = \sigma_g^2

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
                    boost::algorithm::split(parts, line, boost::is_any_of("="));

                    // remove whitespaces from the ends
                    boost::algorithm::trim(parts[0]);

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
};

void Parameter_set::set_value(std::vector<double> vals){
    for (int i=0; i<all.size(); ++i){
        all[i].value =vals[i];
        all[i].miminized = true;
    }
    
}

std::string pad_str(std::string s, const size_t num, const char paddingChar = ' '){
    if(num > s.size())
        s.insert(s.end(), num - s.size(), paddingChar);
    return s;
}

std::ostream& operator<<(std::ostream& os, const Parameter_set& params){
    for (int i=0; i<params.all.size(); ++i){
        if (params.all[i].set){
            if (params.all[i].fixed){
                os <<  pad_str(params.all[i].name, 15) << " (fixed) = " << params.all[i].init;
            } else if (params.all[i].bound){
                os << pad_str(params.all[i].name, 15) << " (bound) = " 
                    << params.all[i].init << ", bounds: (" << params.all[i].lower << ", " << params.all[i].upper << "), step: "
                    << params.all[i].step;
            } else {
                os << pad_str(params.all[i].name, 15) << " (free)  = " 
                    << params.all[i].init << ", step: " << params.all[i].step;
            }
            if (params.all[i].miminized && !params.all[i].fixed){
                os << " -> "<< params.all[i].value;
            }
            os << "\n";
        }
    }
    return os;
}

