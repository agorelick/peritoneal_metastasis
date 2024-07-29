//ad_simulation.cpp
#include <Rcpp.h>
#include <math.h>       /* acos */
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector drift(NumericVector gt, double mu, int gens){
    int n = gt.length();
    int i;
    int j;
    NumericVector variant_markers(n);
    NumericVector indel_direction(n);
    NumericVector gt2(n);

    // copy input genotype to a new vector
    for(j=0; j<n; j++) {
        gt2[j] = gt[j];
    }

    // for each generation, randomly get any indels (presence+direction) and apply them to each marker
    for(i=0; i<gens; i++) {
        variant_markers = Rcpp::rbinom(n, 1, mu);
        indel_direction = Rcpp::rbinom(n, 1, 0.5);
        for(j=0; j<n; j++) {
            if(variant_markers[j]==1 & indel_direction[j]==1) {
                gt2[j] = gt2[j] + 1;
            } else if(variant_markers[j]==1 & indel_direction[j]==0) {
                gt2[j] = gt2[j] - 1;
            } 
        }
    }
    return(gt2);
}

// [[Rcpp::export]]
NumericMatrix angular_distance(NumericMatrix d, int return_Z=0) {
    // get angular distance matrix
    // d should be a matrix of mean-lengths where rows are samples and columns are markers
    // 1st row of d is for the normal

    int nr = d.rows();
    int nc = d.cols();
    int i;
    int j; 
    int s1;
    int s2;
    NumericMatrix Z(nr,nc);
    NumericMatrix AD(nr,nr);

    // subtract normal from every sample
    NumericVector normal = d(0,_);
    for(i=0; i<nr; i++) { 
        d(i,_) = d(i,_) - normal;
    }

    // get Z matrix
    NumericVector denom(nr);
    for(i=0; i<nr; i++) { 
        denom[i] = 0;
        for(j=0; j<nc; j++) { 
            denom[i] = denom[i] + pow(d(i,j),2.0);
        }
        denom[i] = pow(denom[i], 0.5); // sqrt  
        for(j=0; j<nc; j++) {
            Z(i,j) = d(i,j)/denom[i];
        }
    }

    if(return_Z==1) {
        for(j=0; j<nc; j++) {
            Z(0,j) = 0;
        }
        return(Z); 
    }

    // initialize AD matrix with all 0s
    for(s1=0; s1<nr; s1++) { 
        for(s2=0; s2<nr; s2++) { 
            AD(s1,s2) = 0;
        }
    }

    // for each pair of samples in Z, calculate angular distance (exclude excluded samples)
    for(s1=1; s1<nr; s1++) { 
        for(s2=1; s2<nr; s2++) { 
            for(j=0; j<nc; j++) { 
                AD(s1,s2) = AD(s1,s2) + Z(s1,j)*Z(s2,j);
            }
            if(AD(s1,s2) < -1.0) {
                AD(s1,s2) = -1.0;
            }
            if(AD(s1,s2) > 1.0) {
                AD(s1,s2) = 1.0;
            }
            AD(s1,s2) = acos(AD(s1,s2));
            //AD(s2,s1) = AD(s1,s2);
        }
    } 

    // add distances to normal = pi/3
    for(s1=0; s1<nr; s1++) { 
        for(s2=0; s2<nr; s2++) { 
            if((s1==0 | s2==0) & s1!=s2) { 
                AD(s1,s2) = 3.141593/3;
            }
        }
    }

    return(AD);
}
