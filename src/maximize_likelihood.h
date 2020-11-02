

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"
#include <iostream>

#include <math.h> 

double dummt_func(const double *xx ){
    return pow(xx[0] + xx[1], 2);
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
    min->SetPrintLevel(1);

    // create funciton wrapper for minmizer
    // a IMultiGenFunction type
    ROOT::Math::Functor f(&dummt_func,2);
    min->SetFunction(f);


    if (params.mean_lambda.fixed){
        min->SetFixedVariable(0, "mean_lambda", params.mean_lambda.value);
    } else {
        min->SetLimitedVariable(0,"mean_lambda",
                                    params.mean_lambda.value,
                                    params.mean_lambda.step,
                                    params.mean_lambda.lower,
                                    params.mean_lambda.upper);
    }

    if (params.gamma_lambda.fixed){
        min->SetFixedVariable(1, "gamma_lambda",params.gamma_lambda.value);
    } else {
        min->SetLimitedVariable(1,"gamma_lambda",
                                    params.mean_lambda.value,
                                    params.mean_lambda.step,
                                    params.mean_lambda.lower,
                                    params.mean_lambda.upper);
    }


    // do the minimization
    min->Minimize();
    const double *xs = min->X();
    std::cout << "Minimum: f(" << xs[0] << "," << xs[1] << "): " << min->MinValue()  << std::endl;
    std::cout << "Minimizer " << minName << " - " << algoName << std::endl;
}