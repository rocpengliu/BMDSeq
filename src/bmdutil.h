#ifndef BMDUTIL_H
#define BMDUTIL_H

#include <stdlib.h>
#include <string>
#include <iostream>
#include <vector>
#include <sys/stat.h>
#include <algorithm>
#include <time.h>
#include <mutex>
#include <regex>
#include <cmath>
#include <numeric>
#include <iomanip> 
#include "common.h"
#include "util.h"

using namespace std;

struct LinModRes {
    double intercept;
    double slope;
    std::vector<double> fitted_values;
    std::vector<double> residuals;
    double RSS;               // Residual Sum of Squares
    double TSS;               // Total Sum of Squares
    double R_squared;         // Coefficient of Determination
    double std_error_slope;
    double std_error_intercept;
    double AIC;
    size_t n;                 // Number of observations
};

template<typename T>
double get_mean(const std::vector<T>& vec, size_t index, size_t len){
    if (index >= vec.size() || len == 0 || index + len > vec.size()) {
        error_exit("Invalid start or len");
    }
    auto begin = vec.begin() + index;
    auto end = begin + len;
    double sum = std::accumulate(begin, end, 0.0);
    return sum / len;
}

//linear regression
template<typename T1, typename T2>
LinModRes lm(const std::vector<T1>& x, const std::vector<T2>& y, double kcrit = 2.0, int digits = 2){
    size_t n = x.size();
    if(n != y.size() || n == 0){
        error_exit("vectors must be same size and non-empty!");
    }
    double mean_x = get_mean(x, 0, x.size());
    double mean_y = get_mean(y, 0, y.size());
    double Sxy = 0.0;
    double Sxx = 0.0;
    for(size_t i = 0; i < n; ++i){
        Sxy += (x[i] - mean_x)*(y[i] - mean_y);
        Sxx += (x[i] - mean_x) * (x[i] - mean_x);
    }
    double slope = Sxy / Sxx;
    double intercept = mean_y - slope * mean_x;
    std::vector<double> fitted(n), residuals(n);
    double RSS = 0.0, TSS = 0.0;
    for(size_t i = 0; i != n; ++i){
        fitted[i] = intercept + slope * x[i];
        residuals[i] = y[i] - fitted[i];
        RSS += residuals[i] * residuals[i];
        TSS += (y[i] - mean_y) * (y[i] - mean_y);
    }
    double R_squared = 1.0 - (RSS / TSS);
    double variance = RSS / (n - 2);
    double std_error_slope = std::sqrt(variance / Sxx);
    double std_error_intercept = std::sqrt(variance * (1.0 / n + std::pow(mean_x, 2) / Sxx));
    double AIC = kcrit * 2 + n * std::log(RSS / n);

    return {
        intercept,
        slope,
        fitted,
        residuals,
        RSS,
        TSS,
        R_squared,
        std_error_slope,
        std_error_intercept,
        AIC,
        n
    };
}

// template<typename T1, typename T2>
// void fit_model_lin(const std::vector<T1> & x_vec, const std::vector<T2>& dose_vec, double kcrit = 2.0, int digits = 2){
//     auto lin_mod_res = lm(y_vec, dose_vec, kcrit, digits);
//     AIC.i <- round(AIC(fit, k = kcrit), digits = aic_digits)
// }


template<typename T1, typename T2, typename T3>
void perform_bmd(const std::vector<T1> & row_vec, const std::vector<T2>& dose_vec, 
    const std::vector<T3>& dose_rank_vec, double kcrit = 2.0, int digits = 2){
    std::vector<double> y_vec(row_vec.size());
    std::transform(row_vec.begin(), row_vec.end(), y_vec.begin(), [](double x){return std::log2(x+1);});
    auto intercept_slope = lm(y_vec, dose_rank_vec);
    bool adv_incr = intercept_slope.second >= 0;
    for(int i = 0; i != 10; ++i){
        switch(BMD_MODELS[i]){
            case "Lin":
                fit_model_lin(y_vec, dose_vec, kcrit, digits);
                break;

            default:
                break;
        }
    }
}

#endif