functions {

  real single_missing_imputation(real a, real beta, real mu, real sigma){
      // Missing data imputation when one piece of data is missing
      return(Phi((a + beta*mu) / (sqrt(1 + (beta*sigma)^2))));

  }

  row_vector build_observation_process(int state, int y, int ns, real logit_p_detect_int, real logit_p_detect_inf, real rho, int num_states){
    // state: Integer specifying the observed state. Between 1 and 3
    // y: The number of times a frog was detected over the secondary periods
    // ns:  The number of secondary periods
    // logit_p_detect_int:  The log odds of detecting an uninfected frog
    // logit_p_detect_inf:  The change in log odds of detecting an infected frog
    // rho: The probability of a false negative (not load dependent)
    // num_states: For this model, it will be 3 as there are three possible states: Uninfected, Infected, and Dead

    // Detection probability given number of secondary period per primary period
    real dprob_uninf;
    real dprob_inf;
    real p_detect_uninf;
    real p_detect_inf; 
    real tprob;

    // Account for different detection probability between infected and uninfected
    p_detect_uninf = 1 / (1 + exp(-logit_p_detect_int));
    p_detect_inf = 1 / (1 + exp(-(logit_p_detect_int + logit_p_detect_inf)));

    // Account for multiple surveys within primary period
    dprob_uninf = (1 - pow((1 - p_detect_uninf), ns)); // Probability you detect at least once over n surveys
    dprob_inf = (1 - pow((1 - p_detect_inf), ns)); // Probability you detect at least once over n surveys

    // First set of probabilities when true state is Bd -
    // These are probabilities of observing Bd -, Bd + , Unobserved GIVEN true state of Bd-
    if(state == 1){
      tprob = dprob_uninf;  // Given true state of Bd -, probability of observing Bd- is given by detection probability
    } else if (state == 2){
      tprob = 0; // Given true state of Bd-, probability of observing Bd+.  Assuming no False positives
    } else{
      tprob = 1 - dprob_uninf; // Given true state of Bd-, prob of observing undetected
    }
    real pdiag1 = tprob;

    // Second set of probabilities when true state is Bd +
    // These are probabilities of observing Bd -, Bd + , Unobserved GIVEN true state of Bd+
    if(state == 1){
      tprob = dprob_inf*rho;  // Given true state of Bd -, probability of observing Bd- is given by detection probability
    } else if (state == 2){
      tprob = dprob_inf*(1 - rho); // Given true state of Bd-, probability of observing Bd+.  Assuming no False positives
    } else{
      tprob = 1 - dprob_inf; // Given true state of Bd-, prob of observing undetected
    }
    real pdiag2 = tprob;


    // Third set of probabilities when true state is Dead
    // These are the probabilities of observing Bd-, Bd+, unobserved GIVEN true state of Dead
    if(state < num_states){
      tprob = 0; // If you are dead, you can't be observed Bd- or Bd+
    } else{
      tprob = 1;
    }
    real pdiag3 = tprob;

    // Make the diagonal vector
    row_vector[num_states] diag_vector;
    diag_vector[1] = pdiag1;
    diag_vector[2] = pdiag2;
    diag_vector[3] = pdiag3;

    return(diag_vector);

  }



}
data {
  
  int<lower=1> T; // The total number of time steps
  int<lower=1> N; // Number of individuals
  array[N, T] int obs_states; // Observed states where -1 means no individual was available for observation
  array[N, T] int secondary_periods; // The number of secondary periods per primary period
  array[N, T] int secondary_obs; // The number of observations over the secondary periods (always less than or equal to secondary_periods)
  array[N, T] real load; // Loads on individuals
  array[N] int exposure; // Whether individual was exposed or not. 0 unexposed, 1 exposed
  int num_states; // Number of states, 3 for this model
  array[N] int cohort; // 1: for cohort 1 and 2: for cohort 2
  int<lower=1> t2; // The time step at which the second cohort was introduced

}

parameters {

  // Growth parameters
  real mean_load;
  real<lower=0> sigma_load;
  real mean_load_exposure;
  real sigma_load_exposure;
  
  // Loss of infection parameters. Probit
  real loss_int;
  real loss_slope; // Depends on load
  real loss_int_exposure; // Intercept effect of exposure
  
  // Survival parameters. Probit
  real surv_int;
  real<lower=0> surv_slope; // There is not a positive effect of load on survival.
  real surv_int_exposure; // Intercept effect of exposure
  
  // Transmission 
  real trans_int;

  // Detection probability
  real logit_p_detect_int;
  real logit_p_detect_inf;

} 
transformed parameters {

}
model {

  mean_load ~ normal(1, 3); // Allow for large variability in log load
  sigma_load ~ normal(0, 2); // Allow for large variability in sigma

  // Effects of exposure on mean load.  
  mean_load_exposure ~ normal(0, 1.5);
  sigma_load_exposure ~ normal(0, 1.5);

  // On the logit scale, this is weakly centered around 0.5.
  loss_int ~ normal(0, 1.5);

  loss_slope ~ normal(0, 2);

  // Effect of exposure on loss of infection
  loss_int_exposure ~ normal(0, 1);

  // This prior says mean survival probability for uninfected is 0.8 per 7 day
  // time step with mode approaching 1.
  surv_int ~ normal(2, 1.5);
  surv_slope ~ normal(0, 2);

  // Effect of exposure on survival
  surv_int_exposure ~ normal(0, 1.5);

  // On the logit scale, this is weakly centered around 0.5.
  trans_int ~ normal(0, 1.5);

  // Weakly centered about 0.5
  logit_p_detect_int ~ normal(0, 1.5);
  logit_p_detect_inf ~ normal(0, 1.5);

  // Holds shifting quantities
  row_vector[num_states] initial_state; // All individuals start as uninfected upon release
  row_vector[num_states] state_next; // Hold state next
  row_vector[num_states] diag_vector; 
  matrix[num_states, num_states] Tmat = rep_matrix(0, num_states, num_states); // Transition matrix

  // Temporary variables
  real loss_prob;
  real surv_prob;
  real new_surv_int;
  real new_loss_int;
  real new_mean_load;
  real new_sigma_load;
  real marg_loss_prob;
  real marg_surv_prob;


  // Dead stay dead
  Tmat[3, 1] = 0;
  Tmat[3, 2] = 0;
  Tmat[3, 3] = 1;


  // Start as Uninfected.  All individuals start as uninfected
  for(i in 1:num_states){
    if(i == 1){
      initial_state[i] = 1;
    } else{
      initial_state[i] = 0;
    }
  }

  // Loop over individuals
  for(n in 1:N){

    // Allow the effects of loss, survival, and mean to vary with exposure
    if(exposure[n] == 1){
      new_surv_int = surv_int + surv_int_exposure;
      new_loss_int = loss_int + loss_int_exposure;
      new_mean_load = mean_load + mean_load_exposure;
      new_sigma_load = sigma_load*exp(sigma_load_exposure);

    } else{

      new_surv_int = surv_int;
      new_loss_int = loss_int;
      new_mean_load = mean_load;
      new_sigma_load = sigma_load;

    }

    // 
    Tmat[1, 1] = (1 - Phi(trans_int)) * Phi(new_surv_int);
    Tmat[1, 2] = (Phi(trans_int)) * Phi(new_surv_int);
    Tmat[1, 3] = 1 - Phi(new_surv_int);

    // Get the expected loss and survival probability.
    marg_loss_prob = single_missing_imputation(new_loss_int, loss_slope, new_mean_load, new_sigma_load);
    marg_surv_prob = single_missing_imputation(new_surv_int, -1*surv_slope, new_mean_load, new_sigma_load);

    // Loop over an observed trajectory
    for(t in 1:T){

      if(obs_states[n][t] != -1){

        // Survey was conducted
        if (secondary_periods[n][t] != -1){

          // Observation process. Setting \rho = 0. Assuming no observation error as we don't have sufficient power to estimate this.
          diag_vector = build_observation_process(obs_states[n][t], secondary_obs[n][t], secondary_periods[n][t], logit_p_detect_int, logit_p_detect_inf, 0, num_states);

          // Set-up the transition matrix
          if( ((t == 1) && (cohort[n] == 1)) || ((t == t2) && (cohort[n] == 2)) )  {
            // At t = 1 or t = 57, individual is by definition uninfected so the load transitions don't matter. Set to zero
            Tmat[2, 1] = 0;
            Tmat[2, 2] = 0;
            Tmat[2, 3] = 0;

          } else{

            if(obs_states[n][t - 1] == 2){
              // Load was observed on the individual in the previous time step, update the transition probabilities 

              // Update the likelihood for the mean.
              target += normal_lpdf(log(load[n][t - 1]) | new_mean_load, new_sigma_load);

              loss_prob = Phi(new_loss_int + loss_slope*log(load[n][t - 1]));
              surv_prob = Phi(new_surv_int - surv_slope*log(load[n][t - 1] + 1));

              Tmat[2, 1] = loss_prob * surv_prob;
              Tmat[2, 2] = (1 - loss_prob) * surv_prob;
              Tmat[2, 3] = 1 - surv_prob;

            } else{

              // Load was not observed (uninfected or unobserved). Set this to the marginal probs
              Tmat[2, 1] = marg_loss_prob * marg_surv_prob;
              Tmat[2, 2] = (1 - marg_loss_prob) * marg_surv_prob;
              Tmat[2, 3] = (1 - marg_surv_prob);

            }
          }
        
        } else{

          // No survey was conducted
          // Load was not observed (uninfected or unobserved). Set this to the marginal probs
          Tmat[2, 1] = marg_loss_prob * marg_surv_prob;
          Tmat[2, 2] = (1 - marg_loss_prob) * marg_surv_prob;
          Tmat[2, 3] = (1 - marg_surv_prob);

          // We don't need to explicity account for the cohort 2 here because
          // t = 1 from t = 57 will lead to a -1 state and the if statement will
          // never run, which is correct
          if(t != 1){

            if(obs_states[n][t - 1] == 2){

              // Update the likelihood for the mean.
              target += normal_lpdf(log(load[n][t - 1]) | new_mean_load, new_sigma_load);

              loss_prob = Phi(new_loss_int + loss_slope*log(load[n][t - 1]));
              surv_prob = Phi(new_surv_int - surv_slope*log(load[n][t - 1] + 1));

              Tmat[2, 1] = loss_prob * surv_prob;
              Tmat[2, 2] = (1 - loss_prob) * surv_prob;
              Tmat[2, 3] = 1 - surv_prob;
            }

          }
          
          // No survey done
          diag_vector = rep_row_vector(1, num_states);
        }

        // True process x observation process
        if( ((t == 1) && (cohort[n] == 1)) || ((t == t2) && (cohort[n] == 2)) ){
          state_next = (initial_state * Tmat) .* diag_vector;
        } else{
          state_next = (state_next * Tmat) .* diag_vector;
        }

      } // End release period

    } // End T loop

    target += log(sum(state_next));



  } // End individual loop
	
} 

generated quantities {
	
}