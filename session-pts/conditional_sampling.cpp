#include <Rcpp.h>

using namespace Rcpp;

// Function to sample a single conditional survival time given a time grid, survival values, and a censoring time.
// This function uses linear interpolation to find the survival value at the censoring time and then samples
// a time t such that the survival value at t matches the target survival value, which is (1 - uniform_sample) * S(censoring_time).

// [[Rcpp::export]]
List sample_single_conditional_survival_time(
    NumericVector time_grid,
    NumericVector surv_values,
    double censoring_time,
    double tolerance = 1e-6
) {
    int n_times = time_grid.size();
    
    // Linear interpolation for survival value at censoring time.
    double surv_at_cens = 1.0;
    for (int t = 1; t < n_times; t++) {
        if (time_grid[t] >= censoring_time) {
            double t0 = time_grid[t-1];
            double t1 = time_grid[t];
            double s0 = surv_values[t-1];
            double s1 = surv_values[t];
            surv_at_cens = s0 + (s1 - s0) * (censoring_time - t0) / (t1 - t0);
            break;
        }
    }
    
    // Generate a uniform random sample.
    double uniform_sample = Rcpp::runif(1)[0];

    // Target survival value: (1-u) * S(c).
    double target_surv = (1.0 - uniform_sample) * surv_at_cens;
    
    // Find time t where S(t) = target_surv using binary search.
    double t_low = censoring_time;
    double t_high = time_grid[n_times-1];
    double t_result = censoring_time;

    // Be careful with the case where target_surv is smaller than the survival value at the last time point.
    if (target_surv < surv_values[n_times-1]) {
        // This should happen very rarely, but we need to handle it.

        // Let's assume that in this case the result is the last time point.
         return List::create(
            Named("uniform_sample") = uniform_sample,
            Named("t_result") = t_high,
            Named("beyond_max_time") = true
        );
    }
    
    // Binary search for the root.
    for (int iter = 0; iter < 50; iter++) {
        double t_mid = (t_low + t_high) / 2.0;
        
        // Interpolate survival value at t_mid.
        double surv_at_t = 0.0;
        if (t_mid >= time_grid[n_times-1]) {
            surv_at_t = surv_values[n_times-1];
        } else {
            for (int t = 1; t < n_times; t++) {
                if (time_grid[t] >= t_mid) {
                    double t0 = time_grid[t-1];
                    double t1 = time_grid[t];
                    double s0 = surv_values[t-1];
                    double s1 = surv_values[t];
                    surv_at_t = s0 + (s1 - s0) * (t_mid - t0) / (t1 - t0);
                    break;
                }
            }
        }
        
        if (std::abs(surv_at_t - target_surv) < tolerance) {
            t_result = t_mid;
            break;
        }
        
        if (surv_at_t > target_surv) {
            t_low = t_mid;
        } else {
            t_high = t_mid;
        }
        
        t_result = t_mid;
    }
    
    return List::create(
        Named("uniform_sample") = uniform_sample,
        Named("t_result") = t_result,
        Named("beyond_max_time") = false
    );
}

// Function to sample conditional survival times for multiple patients/samples (these are not further distinguished here).
// This function iterates over each sample, calls the single patient sampling function, and collects the
// results into a list of uniform samples and corresponding times.
// It assumes that the input survival values are in a matrix where each row corresponds to a sample
// and each column corresponds to a time point in the time grid.

// [[Rcpp::export]]
List sample_conditional_survival_times(
    NumericVector time_grid, // A vector of time points at which survival values are defined.
    NumericMatrix surv_values, // A matrix of survival values where each row corresponds to a sample and each column corresponds to a time point.
    NumericVector censoring_times, // A vector of censoring times for each sample.
    double tolerance = 1e-6 // Tolerance for convergence in the binary search for the time t.
) {
    int n_samples = surv_values.nrow();
    int n_times = time_grid.size();

    // Assert that the number of censoring times matches the number of samples.
    if (censoring_times.size() != n_samples) {
        stop("Number of censoring times must match number of samples.");
    }

    // Assert that the survival values matrix has the same number of columns as the time grid.
    if (surv_values.ncol() != n_times) {
        stop("Number of columns in survival values matrix must match number of time points in time grid.");
    }
    
    NumericVector uniform_samples(n_samples);
    NumericVector t_results(n_samples);
    LogicalVector beyond_max_time(n_samples);
    
    for (int i = 0; i < n_samples; i++) {
        List result = sample_single_conditional_survival_time(
            time_grid, surv_values.row(i), censoring_times[i], tolerance
        );
        uniform_samples[i] = result["uniform_sample"];
        t_results[i] = result["t_result"];
        beyond_max_time[i] = result["beyond_max_time"];
    }
    
    return List::create(
        Named("uniform_samples") = uniform_samples,
        Named("t_results") = t_results,
        Named("beyond_max_time") = beyond_max_time
    );
}