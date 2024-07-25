//chemo_simulation.cpp
//#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>


using namespace Rcpp;

// [[Rcpp::export]]
NumericVector get_lesion_sizes(NumericMatrix m) { 
    int n_lesions = m.rows();
    int n_clones = m.cols();
    double mysum;
    NumericVector lesion_size(n_lesions);
    
    for(int l=0; l<n_lesions; l++) {
        mysum = 0;
        for(int c=0; c<n_clones; c++) {
            mysum = mysum + m(l,c);
        }
        lesion_size[l] = mysum;
    }
    return(lesion_size); 
}


// [[Rcpp::export]]
NumericMatrix update_population(NumericMatrix m, double prob_death, double max_cells){
    // given a matrix of rows=lesions, columns=clones, where values are number of cells, 
    // update the matrix by one generation, i.e. cells all have a prob of dying (prob_death), and those that don't die all divide.
    // if a lesion has 0 cells of a given clone (extinct) then no longer update it.

    int n_lesions = m.rows();
    int n_clones = m.cols();
    double killed;
    NumericVector lesion_sizes = get_lesion_sizes(m);
 
    for(int l=0; l<n_lesions; l++) {
        if(lesion_sizes[l] > 0 & lesion_sizes[l] < max_cells) {
            for(int c=0; c<n_clones; c++) {
                if(m(l,c) > 0) {
                    killed = R::rbinom(m(l,c), prob_death);
                    m(l,c) = m(l,c) - killed;
                    m(l,c) = 2*m(l,c);
                }
            }
        }
    }

    return(m);
}


// [[Rcpp::export]]
NumericMatrix regrow_lesions(NumericMatrix m, double prob_death, double max_cells, double max_gens=1e5) {

    int l; 
    int c;
    int n_lesions = m.rows();
    int n_clones = m.cols();
    int n_growing = n_lesions;
    IntegerVector growing(n_lesions);
    NumericVector lesion_sizes(n_lesions);

    // initialize a new matrix for the result
    NumericMatrix m2(n_lesions, n_clones);
    for(l=0; l<n_lesions; l++) {
        for(c=0; c<n_clones; c++) {
            m2(l,c) = m(l,c);
        }
    }

    int g=0;
    while(g < max_gens & n_growing > 0) {

        // get the next generation for the entire population
        m2 = update_population(m2, prob_death, max_cells);
        lesion_sizes = get_lesion_sizes(m);

        // check if the lesion is done growing (extinct or max size), if not, update it.
        for(l=0; l<n_lesions; l++) {
            if(lesion_sizes[l] <= 0 | lesion_sizes[l] >= max_cells) {
                growing[l] = 0;
            } else {
                growing[l] = 1;
            }
        }
        n_growing = sum(growing);
        g++;
    }

    return(m2);
}


