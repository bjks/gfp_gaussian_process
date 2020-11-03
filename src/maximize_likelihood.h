

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"
#include <iostream>

#include <math.h> 

double dummt_func(const double *xx ){
    return pow(xx[0] + 2*xx[1], 2);
}

void numerical_minimization(const Parameter_set &params,
                            const char * minName = "Minuit2",
                            const char * algoName = "" ){
    // create minimizer giving a name and a name (optionally) for the specific algorithm

    ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);

    // set tolerance , etc...
    min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
    min->SetMaxIterations(10000);  // for GSL
    min->SetTolerance(0.001);
    min->SetPrintLevel(3);

    // create funciton wrapper for minmizer
    // a IMultiGenFunction type
    ROOT::Math::Functor f(&dummt_func,params.all.size());
    min->SetFunction(f);


    for (int i=0; i<params.all.size(); ++i){
        if (params.all[i].fixed){
            min->SetFixedVariable(i, params.all[i].name, params.all[i].value);
        } else {
            min->SetLimitedVariable(i,  params.all[i].name,
                                        params.all[i].value,
                                        params.all[i].step,
                                        params.all[i].lower,
                                        params.all[i].upper);
        }
    }

 

    // do the minimization
    min->Minimize();
    const double *xs = min->X();
    std::cout << "Minimum: f(" << xs[0] << "," << xs[1] << "): " << min->MinValue()  << std::endl;
    std::cout << "Minimizer " << minName << " - " << algoName << std::endl;
}