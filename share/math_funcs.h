#ifndef math_funcs_h
#define math_funcs_h
#include <math.h>

int logspace_array(unsigned n, double min, double max, double lsa[]){
  unsigned i;
  if ( (max > min) && (min > 0.)){
    double step = (log10(max) - log10(min))/(n-1);
    for (i = 0; i < n; i++){
      lsa[i] = pow(10.,log10(min) + i*step);
      }
    return 0;
    } else {
    return 1;
    }
  }

int linspace_array(unsigned n, double min, double max, double lsa[]){
  unsigned i;
  if ( (max > min) && (min > 0.)){
    double step = (max - min)/(n-1);
    for (i = 0; i < n; i++){
      lsa[i] = min + i*step;
      }
    return 0;
    } else {
    return 1;
    }
  }

#endif
