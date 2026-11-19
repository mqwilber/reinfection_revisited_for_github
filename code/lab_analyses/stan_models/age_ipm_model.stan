data{
    // int D; // Number of observations
    int<lower=0> T; // Number of treatments
    int<lower=0> H; // Number of hosts per treatment
    int<lower=0> time_steps; // Number of time steps
    real samps[T, H, time_steps]; // 3d matrix to hold swabbing time series

} parameters{

    // Growth function parameters
    real bg0; // intercept
    real bg1; // fungal effect
    real bg_tA; // Effect of recent infected relative to naive on growth on slope
    real bg_tA2; // Effect of recent infected relative to naive on growth on intercept


    // Initial infection function parameters
    real binit0; // intercept

    // Variance parameters
    real<lower=0> sigma_init; // Standard deviation of initial infection?
    real<lower=0> sigma_growth; // Variance of growth


} model {
    
    real prob_inf; // Probability of infection in a time step
    real mean_init_inf; // Mean initial infection in a time step
    real mean_growth; // Mean Bd growth
    real loss_prob; // Loss of Bd infection
    int treat_coef; // Holds treatment coefs
    real power_sd;

    // Priors on parameters. 
    bg1 ~ normal(0, 5);
    bg_tA ~ normal(0, 1);
    bg_tA2 ~ normal(0, 1); 

    sigma_init ~ cauchy(0, 2);
    sigma_growth ~ cauchy(0, 2);

    // Loop over treatments
    for(t in 1:T) { 

        // Assign treatment effects
        if(t == 1){
            treat_coef = 0;
        } else{
            treat_coef = 1;
        } 

        // Loop over hosts per treatment
        for(h in 1:H) {

            // Loop over the number of time steps
            for(time in 1:time_steps) { 

                // Only consider data points that exist, because we can not use ragged arrays

                if((samps[t, h, time] != -1) && (samps[t, h, time + 1] != -1)) { 

                    // Choose the right likelihood increment

                    // 1. Transmission and initial infection
                    if(samps[t, h, time] == 0 && samps[t, h, time + 1] > 0){
                         

                        mean_init_inf = binit0; // - binit1*i_current; 

                        target += normal_lpdf(log(samps[t, h, time + 1]) | mean_init_inf, sigma_init);
                    }

                    // 2. Dont lose infection and growth of infection
                    if((samps[t, h, time] > 0) && (samps[t, h, time + 1] > 0)){
                        
                        mean_growth = bg0 + bg_tA2*treat_coef + log(samps[t, h, time])*(bg1 + bg_tA*treat_coef); 
                        target += normal_lpdf(log(samps[t, h, time + 1]) | mean_growth, sigma_growth);

                    }

                }

            }

        }

    }   

}   