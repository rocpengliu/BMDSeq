#ifndef BMD_SEEKER_H
#define BMD_SEEKER_H

#include "bmdscore/bmds_helper.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

enum model_type { dichotomous, continuous};

using namespace std;

class BmdSeeker{
public:
    BmdSeeker();
    ~BmdSeeker();

public:
    void runOldDichoAnalysis();
    void runOldContAnalysis();//right one
    void runDichoMA();
    void runPythonDichoAnalysis();
    void runPythonDichoMA();
    void runPythonContAnalysis(const std::vector<double> & x_vec, 
                               const std::vector<double> & y_vec, bool trend);//right one
    void runPythonMultitumorAnalysis();
    void runTestMultitumorModel();
    void runPythonNestedAnalysis();
    void test();
    void printDichoModResult();
    void printDichoModResult(struct python_dichotomous_analysis *pyAnal, 
                                        struct python_dichotomous_model_result *pyRes, 
                                        bool showResultsOverride);
    void printNestedModResult(struct python_nested_analysis *pyAnal, 
                              struct python_nested_result *pyRes,
                              bool showResultsOverride);
    void printNestedModResult();
    void Nlogist_probs_test();
    void Nlogist_lk_test();
    void runTestMTInequalityConstraint();
    void runTestMTEqualityConstraint();

public:
    std::vector<double> getMultitumorPrior(int degree, int prior_cols);
    std::vector<double> getNlogisticPrior(int ngrp, int prior_cols, bool restricted);
    std::vector<double> getNctrPrior(int ngrp, int prior_cols, bool restricted);
    bool showResultsOverride;

};

#endif