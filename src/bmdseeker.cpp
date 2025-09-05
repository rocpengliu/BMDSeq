#include "bmdseeker.h"
#include "priors.h"
#include <chrono>

BmdSeeker::BmdSeeker(){
  showResultsOverride = true;
}

BmdSeeker::~BmdSeeker(){

}

void BmdSeeker::runPythonContAnalysis(const std::vector<double> & x_vec, 
                                      const std::vector<double> & y_vec, bool trend) {
  //printf("Running python continuous analysis\n");
  auto start = std::chrono::high_resolution_clock::now();

  bool isIncreasing;

  ///////////////////////////////
  // USER INPUT
  //////////////////////////////

  enum cont_model model = exp_5;    // hill, exp_3, exp_5, power, funl, polynomial
  int modelType = 1;                // 1 = frequentist, 2 = bayesian
  bool restricted = true;           // only used for frequentist models
  enum distribution dist = normal;  // normal, normal_ncv, log_normal
  bool detectAdvDir = true;         // if false then need to set isIncreasing
  bool penalizeAIC = true;
  // isIncreasing = true;

  int degree = 2;  // for polynomial only

  double alpha = 0.05;
  double BMRF = 1.0;  // 1.0;
  int BMD_type = 2;   // 1=absdev, 2 = stddev, 3 = reldev, 4 = pt, 5 = extra, 6 = hybrid_extra, 7 =
                      // hybrid_added   from src/include/cmodeldefs.h
  ////////////////////////////////////////////////
  // cont data - suff stat: dose, Y, N, SD
  // cont data - individual: dose, response
  /////////////////////////////////////////////
  bool suffStat = false;

  // continuous2.dax
  // double* D = new double[x_vec.size()];
  // double* Y = new double[y_vec.size()];
  // for(int i = 0; i != x_vec.size(); ++i){
  //   D[i] = x_vec[i];
  //   Y[i] = static_cast<double>(y_vec[i]);
  // }
  // for(int i = 0; i != x_vec.size(); ++i){
  //   std::cout << "index: " << i << " : " << D[i] << " DOSE " << Y[i] << " response\n";
  // }
  // double D[] = {0,  0,  0,  0,  18, 18, 18, 18, 18, 20, 20, 20, 20,
  //               30, 30, 30, 30, 35, 35, 35, 35, 40, 40, 40, 40, 40};
  // double Y[] = {39, 38.4, 36.3, 37.1, 40.2, 45.3, 42.1, 38.3, 35.9, 42.5, 45.2, 40.1, 39.8,
  //               50.1, 53.4, 48.2, 52.1, 56.1, 50.4, 53.2, 55.2, 55.1, 59.1, 56.3, 52.9, 53.7};
  double D[] = {0,  0,  0, 
                0.000001, 0.000001, 0.000001,
                0.00001, 0.00001, 0.00001,
                0.0003, 0.0003, 0.0003,
                0.001, 0.001, 0.001,
                0.003, 0.003, 0.003,
                0.01, 0.01, 0.01,
                0.03, 0.03, 0.03,
                0.1, 0.1, 0.1,
                1, 1, 1};
  double Y[] = {10131, 10812, 9413,
                9693, 11100, 10171,
                8974, 10944, 10065,
                11051, 10104, 9507, 
                9668, 10692, 7883, 
                9810, 10532, 8211,
                9945, 9653, 8859, 
                9907, 9256, 9892,
                10790, 9757, 10240,
                8934, 12281, 10716};
  double N[1];
  double SD[1];
  //isIncreasing = trend;
  isIncreasing = true;

  /////////////////////////////////////////////////
  // END USER INPUT
  ///////////////////////////////////////////////////

  // struct continuous_analysis anal;
  struct python_continuous_analysis anal;
  anal.penalizeAIC = penalizeAIC;
  int numDataRows = sizeof(D) / sizeof(D[0]);
  //int numDataRows = x_vec.size();

  if (!detectAdvDir) {
    anal.isIncreasing = isIncreasing;
  }

  // check data array sizes for consistency
  size_t numElementsY = sizeof(Y) / sizeof(Y[0]);
  //size_t numElementsY = y_vec.size();

  if (!detectAdvDir) {
    anal.isIncreasing = isIncreasing;
  }
  if (suffStat) {
    size_t numElementsN = sizeof(N) / sizeof(N[0]);
    size_t numElementsSD = sizeof(SD) / sizeof(SD[0]);
    if (numDataRows != numElementsY || numElementsY != numElementsN ||
        numElementsN != numElementsSD) {
      printf("Number of data elements are not consistent\nExiting Code\n");
      exit(-1);
    }
  } else {
    if (numDataRows != numElementsY) {
      printf("Number of data elements are not consistent\nExiting Code\n");
      exit(-1);
    }
  }

  // priors defined columnwise
  int prCols = 5;

  // define priors/parameter constraints
  int numParms;
  //printf("model = %d\n", model);
  switch (model) {
    case hill:
      numParms = 6;
      break;
    case exp_3:
      // numParms = 5;
      // break;
    case exp_5:
      numParms = 6;
      break;
    case power:
      numParms = 5;
      break;
    case funl:
      numParms = 0;  // FIX THIS
      break;
    case polynomial:
      numParms = 3 + degree;
      break;
    default:
      printf("error in numParms\n");
      return;
  }
  if (dist == normal || dist == log_normal) {
    numParms -= 1;
  }

  //printf("numParms = %d\n", numParms);
  // double* pr;
  double *prior;

  //printf("starting priors\n");

  if (modelType == 1) {
    // frequentist
    if (restricted) {
      printf("choosing frequentist restricted priors\n");
      switch (model) {
        case hill:
          anal.model = hill;
          if (dist == normal || dist == log_normal) {
            // normal
            prior = prRFreqHillNormal;
          } else {
            //} else if (dist == normal_ncv){
            // normal NCV
            prior = prRFreqHillNormalNCV;
            //} else {
            //  //lognormal
            //  anal.prior = prRFreqHillLognormal;
          }
          break;
        case exp_3:
          anal.model = exp_3;
          //          if (dist == normal || dist == log_normal){
          //            anal.prior = prRFreqExp5Normal;
          //          } else {
          //            anal.prior = prRFreqExp5NormalNCV;
          //          }
          if (dist == normal) {
            prior = prRFreqExp5Normal;
          } else if (dist == normal_ncv) {
            prior = prRFreqExp5NormalNCV;
          } else {
            prior = prRFreqExp5Lognormal;
          }
          break;
        case exp_5:
          anal.model = exp_5;
          //          if (dist == normal || dist == log_normal){
          //            anal.prior = prRFreqExp5Normal;
          //          } else {
          //            anal.prior = prRFreqExp5NormalNCV;
          //          }
          if (dist == normal) {
            prior = prRFreqExp5Normal;
          } else if (dist == normal_ncv) {
            prior = prRFreqExp5NormalNCV;
          } else {
            prior = prRFreqExp5Lognormal;
          }
          break;
        case power:
          anal.model = power;
          if (dist == normal || dist == log_normal) {
            prior = prRFreqPower;
          } else {
            prior = prRFreqPowerNCV;
          }
          break;
        case funl:
          break;
        case polynomial:
          printf("choosing polynomial model\n");
          anal.model = polynomial;
          anal.degree = degree;
          if (detectAdvDir) {
            if (dist == normal || dist == log_normal) {
              printf("using advDir auto normal or log_normal dist priors\n");
              if (degree == 1) {
                prior = prRFreqPoly1;
              } else if (degree == 2) {
                prior = prRFreqPoly2;
              } else if (degree == 3) {
                prior = prRFreqPoly3;
              } else if (degree == 4) {
                prior = prRFreqPoly4;
              } else if (degree == 5) {
                prior = prRFreqPoly5;
              } else {
                printf("poly restricted normal/lognormal degree error\n");
                return;
              }
            } else {
              printf("using advDir auto normal_ncv dist priors\n");
              if (degree == 1) {
                prior = prRFreqPoly1NCV;
              } else if (degree == 2) {
                prior = prRFreqPoly2NCV;
              } else if (degree == 3) {
                prior = prRFreqPoly3NCV;
              } else if (degree == 4) {
                prior = prRFreqPoly4NCV;
              } else if (degree == 5) {
                prior = prRFreqPoly5NCV;
              } else {
                printf("poly restricted normal NCV degree error\n");
                return;
              }
            }
          } else {
            if (anal.isIncreasing) {
              if (dist == normal || dist == log_normal) {
                printf("using advDir up normal or log_normal dist priors\n");
                if (degree == 1) {
                  prior = prRFreqPoly1Up;
                } else if (degree == 2) {
                  prior = prRFreqPoly2Up;
                } else if (degree == 3) {
                  prior = prRFreqPoly3Up;
                } else if (degree == 4) {
                  prior = prRFreqPoly4Up;
                } else if (degree == 5) {
                  prior = prRFreqPoly5Up;
                } else {
                  printf("poly restricted normal/lognormal degree error\n");
                  return;
                }
              } else {
                printf("using advDir up normal_ncv dist priors\n");
                if (degree == 1) {
                  prior = prRFreqPoly1NCVUp;
                } else if (degree == 2) {
                  prior = prRFreqPoly2NCVUp;
                } else if (degree == 3) {
                  prior = prRFreqPoly3NCVUp;
                } else if (degree == 4) {
                  prior = prRFreqPoly4NCVUp;
                } else if (degree == 5) {
                  prior = prRFreqPoly5NCVUp;
                } else {
                  printf("poly restricted normal NCV degree error\n");
                  return;
                }
              }

            } else {
              if (dist == normal || dist == log_normal) {
                printf("using advDir down normal or log_normal dist priors\n");
                if (degree == 1) {
                  prior = prRFreqPoly1Down;
                } else if (degree == 2) {
                  prior = prRFreqPoly2Down;
                } else if (degree == 3) {
                  printf("using prRFreqPoly3Down\n");
                  prior = prRFreqPoly3Down;
                } else if (degree == 4) {
                  prior = prRFreqPoly4Down;
                } else if (degree == 5) {
                  prior = prRFreqPoly5Down;
                } else {
                  printf("poly restricted normal/lognormal degree error\n");
                  return;
                }
              } else {
                printf("using advDir down normal_ncv dist priors\n");
                if (degree == 1) {
                  prior = prRFreqPoly1NCVDown;
                } else if (degree == 2) {
                  prior = prRFreqPoly2NCVDown;
                } else if (degree == 3) {
                  prior = prRFreqPoly3NCVDown;
                } else if (degree == 4) {
                  prior = prRFreqPoly4NCVDown;
                } else if (degree == 5) {
                  prior = prRFreqPoly5NCVDown;
                } else {
                  printf("poly restricted normal NCV degree error\n");
                  return;
                }
              }
            }
          }
          break;
        default:
          printf("error with restricted models\n");
          return;
      }
    } else {
      // unrestricted
      switch (model) {
        case hill:
          anal.model = hill;
          if (dist == normal) {
            // normal
            prior = prUFreqHillNormal;
          } else if (dist == normal_ncv) {
            // normal NCV
            prior = prUFreqHillNormalNCV;
          } else {
            // lognormal
            prior = prUFreqHillLognormal;
          }
          break;
        case exp_3:
          printf("cannot run unrestricted exponential models\n");
          return;
          // break;
        case exp_5:
          printf("cannot run unrestricted exponential models\n");
          return;
          // break;
        case power:
          anal.model = power;
          if (dist == normal || dist == log_normal) {
            prior = prUFreqPower;
          } else {
            prior = prUFreqPowerNCV;
          }
          break;

        case funl:
          break;
        case polynomial:
          printf("choosing polynomial model\n");
          anal.model = polynomial;
          anal.degree = degree;
          // if (detectAdvDir){
          if (dist == normal || dist == log_normal) {
            printf("prior with normal or lognormal dist\n");
            if (degree == 1) {
              prior = prUFreqPoly1;
            } else if (degree == 2) {
              prior = prUFreqPoly2;
            } else if (degree == 3) {
              prior = prUFreqPoly3;
            } else if (degree == 4) {
              prior = prUFreqPoly4;
            } else if (degree == 5) {
              prior = prUFreqPoly5;
            } else {
              printf("poly unrestricted normal/lognormal degree error\n");
              return;
            }
          } else {
            if (degree == 1) {
              prior = prUFreqPoly1NCV;
            } else if (degree == 2) {
              prior = prUFreqPoly2NCV;
            } else if (degree == 3) {
              prior = prUFreqPoly3NCV;
            } else if (degree == 4) {
              prior = prUFreqPoly4NCV;
            } else if (degree == 5) {
              prior = prUFreqPoly5NCV;
            } else {
              printf("poly restricted normal NCV degree error\n");
              return;
            }
          }
          //}
          break;

        default:
          printf("error with unrestricted model\n");
          return;
      }
    }
  } else {
    // bayesian
    switch (model) {
      case hill:
        anal.model = hill;
        if (dist == normal || dist == log_normal) {
          // normal
          prior = prBayesianHill;
        } else {
          // normal NCV
          prior = prBayesianHillNCV;
        }
        break;
      case exp_3:
        anal.model = exp_3;
        if (dist == normal || dist == log_normal) {
          // normal
          prior = prBayesianExp5;
        } else {
          // normal NCV
          prior = prBayesianExp5NCV;
        }
        break;
      case exp_5:
        anal.model = exp_5;
        if (dist == normal || dist == log_normal) {
          // normal
          prior = prBayesianExp5;
        } else {
          // normal NCV
          prior = prBayesianExp5NCV;
        }
        break;
      case power:
        anal.model = power;
        if (dist == normal || dist == log_normal) {
          // normal
          prior = prBayesianPower;
        } else {
          // normal NCV
          prior = prBayesianPowerNCV;
        }
        break;
      case funl:
        anal.model = funl;
        printf("FUNL model has not been implemented in BMDS\n");
        break;
      case polynomial:
        anal.model = polynomial;
        anal.degree = degree;
        if (dist == normal || dist == log_normal) {
          // normal
          printf("using Bayesian normal or log_normal dist priors\n");
          if (degree == 1) {
            prior = prBayesianPoly1;
          } else if (degree == 2) {
            prior = prBayesianPoly2;
          } else if (degree == 3) {
            prior = prBayesianPoly3;
          } else if (degree == 4) {
            prior = prBayesianPoly4;
          } else if (degree == 5) {
            prior = prBayesianPoly5;
          } else {
            printf("poly restricted normal/lognormal degree error\n");
            return;
          }
        } else {
          // normal NCV
          printf("using Bayesian normal_ncv dist priors\n");
          if (degree == 1) {
            prior = prBayesianPoly1NCV;
          } else if (degree == 2) {
            prior = prBayesianPoly2NCV;
          } else if (degree == 3) {
            prior = prBayesianPoly3NCV;
          } else if (degree == 4) {
            prior = prBayesianPoly4NCV;
          } else if (degree == 5) {
            prior = prBayesianPoly5NCV;
          } else {
            printf("poly restricted normal/lognormal degree error\n");
            return;
          }
        }
        break;
    }
  }

  printf("initial priors\n");
  //  for (int i=0; i<numParms * anal.prior_cols; i++){
  //    printf("%f,",anal.prior[i]);
  //  }

  printf("finished with priors\n");
  //

  // parms array declared
  //  int numParms = sizeof(pr)/sizeof(pr[0])/prCols;
  // double parms[numParms];
  double *parms = new double[numParms];

  // declare analysis
  anal.Y.assign(Y, Y + numDataRows);
  //anal.Y.assign(y_vec.begin(), y_vec.end());
  anal.n = numDataRows;
  if (suffStat) {
    anal.n_group.assign(N, N + numDataRows);
    anal.sd.assign(SD, SD + numDataRows);
  }
  anal.doses.assign(D, D + numDataRows);
  //anal.doses.assign(x_vec.begin(), x_vec.end());
  anal.disttype = dist;
  if (!detectAdvDir) {
    anal.isIncreasing = isIncreasing;
  }
 printf("1111111111111111111111111\n");
  anal.alpha = alpha;
  anal.BMD_type = BMD_type;  // 1=absdev, 2 = stddev, 3 = reldev, 4 = pt, 5 = extra, 6 =
                             // hybrid_extra, 7 = hybrid_added   from src/include/cmodeldefs.h
  anal.BMR = BMRF;
  anal.samples = 0;  // num MCMC samples
  anal.tail_prob = 0.01;
  anal.suff_stat = suffStat;
  anal.parms = numParms;
  anal.prior_cols = prCols;
  anal.transform_dose = 0;
  anal.prior.assign(prior, prior + anal.prior_cols * anal.parms);
  anal.restricted = restricted;
  anal.detectAdvDir = detectAdvDir;

  //printf("prior b4 adj:\n");
  // for (int i = 0; i < prCols * numParms; i++) {
  //   printf("%.9f\n", anal.prior[i]);
  // }
 printf("222222222222222222222\n");
  struct python_continuous_model_result res;
  res.model = anal.model;
  res.nparms = anal.parms;
  res.dist_numE = 100;

  struct BMDS_results BMDSres;
  // set all parms as unbounded initially
  for (int i = 0; i < anal.parms; i++) {
    BMDSres.bounded.push_back(false);
    BMDSres.stdErr.push_back(BMDS_MISSING);
    BMDSres.lowerConf.push_back(BMDS_MISSING);
    BMDSres.upperConf.push_back(BMDS_MISSING);
  }
  BMDSres.BMD = -9999.0;
  BMDSres.BMDU = -9999.0;
  BMDSres.BMDL = -9999.0;
  BMDSres.AIC = -9999.0;
  res.bmdsRes = BMDSres;

  struct continuous_GOF gof;
  int nGOF;
  if (anal.suff_stat) {
    nGOF = anal.n;
  } else {
    // determine number of unique dose groups
    nGOF = 1;
    for (int i = 1; i < anal.n; i++) {
      int j = 0;
      for (j = 0; j < i; j++) {
        if (anal.doses[i] == anal.doses[j]) break;
      }
      if (i == j) nGOF++;
    }
  }
  res.gof = gof;
 printf("3333333333333333333333333333\n");
  struct continuous_AOD aod;
  std::vector<double> LL(5);
  std::vector<int> nParms(5);
  std::vector<double> AIC(5);
  double addConst;  // = 22.2;
  aod.LL = LL;
  aod.nParms = nParms;
  aod.AIC = AIC;
  aod.addConst = addConst;

  // double llRatio[4];
  // double DF[4];
  // double pVal[4];
  std::vector<double> llRatio(4);
  std::vector<double> DF(4);
  std::vector<double> pVal(4);
  struct testsOfInterest TOI;

  TOI.llRatio = llRatio;
  TOI.DF = DF;
  TOI.pVal = pVal;
  aod.TOI = TOI;
  res.aod = aod;

  printf("\n\n-----------INPUT---------------\n");
  //printBmdsStruct(&anal);
  //printf("\n\n");

  printf("calling pythonBMDSCont\n");
  pythonBMDSCont(&anal, &res);

  //printf("\n\n");
  printf("prior after adj by model code:\n");
  // for (int i = 0; i < prCols * numParms; i++) {
  //   printf("%.20f\n", anal.prior[i]);
  // }

  printf("\n\n----------OUTPUT-----------\n");
  if (BMDSres.validResult || showResultsOverride) {
    //printBmdsStruct(&res);
    //    printf("\nBMD Dist:\n");
    //    for (int i = 0; i < res.dist_numE; i++) {
    //      printf("i:%d, perc:%f, dist:%f\n", i, res.bmd_dist[i + res.dist_numE], res.bmd_dist[i]);
    //    }
    //
  } else {
    //printf("Model was not run\n");
  }

  // delete[] D;
  // delete[] Y;
  auto end = std::chrono::high_resolution_clock::now();
  double duration = std::chrono::duration<double>(end-start).count();
  printf("Function executed in %.6f seconds\n", duration);
  delete[] parms;
  // debugContAnal(&anal);
}