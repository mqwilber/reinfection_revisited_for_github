data{
    int D; // Number of observations
    int<lower=0> T; // Number of treatments
    int<lower=0> H; // Number of hosts per treatment
    int<lower=0> time_steps; // Number of time steps
    real samps[T, H, time_steps]; // 3d matrix to hold swabbing time series

} parameters{

    // Growth function parameters
    real bg0; // intercept
    real bg1; // fungal effect
    real bg_tA; // Effect of recent infected relative to naive on growth on slope
    real bg_tB; // Effect of distance infected relative to naive on growth on slope
    real bg_tA2; // Effect of recent infected relative to naive on growth on intercept
    real bg_tB2; // Effect of distance infected relative to naive on growth on intercept

    // Loss of infection function parameters
    real bl0; // intercept
    real bl1; // fungal effect
    real bl_tA; // Effect of recently infected relative to naive
    real bl_tB; // Effect of distant infected relative to naive
    // real<lower=0> bl3; // immunity effect squared

    // Transmission function parameters
    real bt0; // intercept for transmission after first dose

    // Transmission function for first dose
    real bt_init0; // intercept for first dose infection (effect on naive)
    real bt_tA; // Effect of recently infected relative to naive 
    real bt_tB; // Effect of distant infected relative to naive

    // real<lower=0> bt1; // immunity effect

    // Initial infection function parameters
    real binit0; // intercept

    // Variance parameters
    real<lower=0> sigma_init; // Standard deviation of initial infection?
    real<lower=0> sigma2_growth; // Variance of growth
    real delta; // Allowing variance in Bd growth to change with with infection load


} model {
    
    real prob_inf; // Probability of infection in a time step
    real mean_init_inf; // Mean initial infection in a time step
    real mean_growth; // Mean Bd growth
    real loss_prob; // Loss of Bd infection
    int treat_coef[2]; // Holds treatment coefs
    real power_sd;

    // Priors on parameters. 
    bg1 ~ normal(0, 5);
    bg_tA ~ normal(0, 1); // Very little information to estimate A
    bg_tB ~ normal(0, 5);
    bg_tA2 ~ normal(0, 1); // Very little information to estimate A
    bg_tB2 ~ normal(0, 5);

    bl1 ~ normal(0, 5); 
    bl_tA ~ normal(0, 5); 
    bl_tB ~ normal(0, 5);

    bt_init0 ~ normal(0, 5);
    bt_tA ~ normal(0, 5);
    bt_tB ~ normal(0, 5);

    sigma_init ~ cauchy(0, 2);
    sigma2_growth ~ cauchy(0, 2);
    delta ~ normal(0, 5);

    // Loop over treatments
    for(t in 1:T) { 

        // Assign treatment effects
        if(t == 1){
            treat_coef[1] = 0;
            treat_coef[2] = 0;
        } else if(t == 2){
            treat_coef[1] = 1;
            treat_coef[2] = 0;
        } else{
            treat_coef[1] = 0;
            treat_coef[2] = 1;
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
                        
                        // Account for infection at initial dose
                        if(time == 1)
                            prob_inf = inv_logit(bt_init0 + bt_tA*treat_coef[1] + bt_tB*treat_coef[2]);
                        else
                            prob_inf = inv_logit(bt0);

                        mean_init_inf = binit0; // - binit1*i_current; 

                        target += bernoulli_lpmf(1 | prob_inf) + normal_lpdf(log(samps[t, h, time + 1]) | mean_init_inf, sigma_init);
                    }

                    // 2. Dont lose infection and growth of infection
                    if((samps[t, h, time] > 0) && (samps[t, h, time + 1] > 0)){
                        
                        loss_prob = inv_logit(bl0 + bl1*log(samps[t, h, time]) + bl_tA*treat_coef[1] + bl_tB*treat_coef[2]);
                        power_sd = sqrt(sigma2_growth * exp(2*delta*log(samps[t, h, time])));
                        mean_growth = bg0 + log(samps[t, h, time])*(bg1 + bg_tA*treat_coef[1] + bg_tB*treat_coef[2]) + bg_tA2*treat_coef[1] + bg_tB2*treat_coef[2]; // CHECK IMMUNITY EFFECT
                        target += bernoulli_lpmf(0 | loss_prob) + normal_lpdf(log(samps[t, h, time + 1]) | mean_growth, power_sd);

                    }

                    // 3. Loss of infection
                    if((samps[t, h, time] > 0) && (samps[t, h, time + 1] == 0)){
                        
                        loss_prob = inv_logit(bl0 + bl1*log(samps[t, h, time]) + bl_tA*treat_coef[1] + bl_tB*treat_coef[2]);
                        target += bernoulli_lpmf(1 | loss_prob);
                    }

                    // 4. No transmission
                    if((samps[t, h, time] == 0) && (samps[t, h, time + 1] == 0)){
                        
                        // Account for initial dose
                        if(time == 1)
                            prob_inf = inv_logit(bt_init0 + bt_tA*treat_coef[1] + bt_tB*treat_coef[2]);
                        else
                            prob_inf = inv_logit(bt0); // - bt1*i_current);

                        target += bernoulli_lpmf(0 | prob_inf);
                    }

                }

            }

        }

    }   

} generated quantities {

    // Extracting the log-likelihoods for model comparison purposes
    vector[D] log_lik;
    int counter;
    real dev;
    real prob_inf; // Probability of infection in a time step
    real mean_init_inf; // Mean initial infection in a time step
    real mean_growth; // Mean Bd growth
    real loss_prob; // Loss of Bd infection
    real power_sd;
    int treat_coef[2];
    // real immunity[T, H, time_steps - 1]; // Hold immunity projections over time

    counter = 0;
    dev = 0;

   // Loop over treatments
    for(t in 1:T) { 

        if(t == 1){
            treat_coef[1] = 0;
            treat_coef[2] = 0;
        } else if(t == 2){
            treat_coef[1] = 1;
            treat_coef[2] = 0;
        } else{
            treat_coef[1] = 0;
            treat_coef[2] = 1;
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
                        
                        // Account for infection at initial dose
                        if(time == 1)
                            prob_inf = inv_logit(bt_init0 + bt_tA*treat_coef[1] + bt_tB*treat_coef[2]);
                        else
                            prob_inf = inv_logit(bt0);

                        mean_init_inf = binit0; // - binit1*i_current; 

                        counter = counter + 1;
                        log_lik[counter] = bernoulli_lpmf(1 | prob_inf) + normal_lpdf(log(samps[t, h, time + 1]) | mean_init_inf, sigma_init);
                        dev = dev - (2 * log_lik[counter]);
                    }

                    // 2. Dont lose infection and growth of infection
                    if((samps[t, h, time] > 0) && (samps[t, h, time + 1] > 0)){

                        loss_prob = inv_logit(bl0 + bl1*log(samps[t, h, time]) + bl_tA*treat_coef[1] + bl_tB*treat_coef[2]);
                        power_sd = sqrt(sigma2_growth * exp(2*delta*log(samps[t, h, time])));
                        mean_growth = bg0 + log(samps[t, h, time])*(bg1 + bg_tA*treat_coef[1] + bg_tB*treat_coef[2]) + bg_tA2*treat_coef[1] + bg_tB2*treat_coef[2];

                        counter = counter + 1;
                        log_lik[counter] = bernoulli_lpmf(0 | loss_prob) + normal_lpdf(log(samps[t, h, time + 1]) | mean_growth, power_sd);
                        dev = dev - (2 * log_lik[counter]);
                    }

                    // 3. Loss of infection
                    if((samps[t, h, time] > 0) && (samps[t, h, time + 1] == 0)){
                        
                        loss_prob = inv_logit(bl0 + bl1*log(samps[t, h, time]) + bl_tA*treat_coef[1] + bl_tB*treat_coef[2]);
                        counter = counter + 1;
                        log_lik[counter] = bernoulli_lpmf(1 | loss_prob);
                        dev = dev - (2 * log_lik[counter]);
                    }

                    // 4. No transmission
                    if((samps[t, h, time] == 0) && (samps[t, h, time + 1] == 0)){
                        
                        // Account for initial dose
                        if(time == 1)
                            prob_inf = inv_logit(bt_init0 + bt_tA*treat_coef[1] + bt_tB*treat_coef[2]);
                        else
                            prob_inf = inv_logit(bt0); // - bt1*i_current);
                        
                        counter = counter + 1;
                        log_lik[counter] = bernoulli_lpmf(0 | prob_inf);
                        dev = dev - (2 * log_lik[counter]);
                    }

                }
            }

        }

    }
}   
   