data {
  // experiment parameters
  real stepSec;// duration between two decision points
  real iti;// iti duration
  int  nWaitOrQuit; // number of possible decision points
  real tWaits[nWaitOrQuit]; // time for each decision point 
  
  // initial value for V0
  real V0_ini; 
  
  // empirical data
  int N; // number of trials
  int Rs[N]; // payoff in each trial
  real Ts[N]; // a trial ends at t == T
  int nMadeActions[N];// at which decision point each trial ends
}
transformed data {
  // total number of decision points in all trials
  int nTotalAction = sum(nMadeActions);
}
parameters {
  // parameters:
  // alphaR : learning rate for rewards
  // alphaU : learning rate for non-rewards
  // tau : action consistency, namely the soft-max temperature parameter
  // gamma: discount factor
  // prior: prior belief parameter
  
  // for computational efficiency,we sample raw parameters from unif(-0.5, 0.5)
  // which are later transformed into actual parameters
  real<lower = -0.5, upper = 0.5> raw_alpha;
  real<lower = -0.5, upper = 0.5> raw_rho; // valence-based bias
  real<lower = -0.5, upper = 0.5> raw_tau;
  real<lower = -0.5, upper = 0.5> raw_gamma;
  real<lower = -0.5, upper = 0.5> raw_eta;
}
transformed parameters{
  // scale raw parameters into real parameters
  real alpha = (raw_alpha + 0.5) * 0.3; // alphaR ~ unif(0, 0.3)
  real alphaU = min([alpha * (raw_rho + 0.5) * 5, 1]');// alphaU
  real rho = alphaU / alpha;
  real tau = (raw_tau + 0.5) * 21.9 + 0.1; // tau ~ unif(0.1, 22)
  real gamma = (raw_gamma + 0.5) * 0.3 + 0.7; // gamma ~ unif(0.7, 1)
  real eta = (raw_eta + 0.5) * 6.5; // prior ~ unif(0, 6.5)
  
  // declare variables 
  // // state value of t = 0
  real V0; 
  // // action value of waiting in each decision points 
  vector[nWaitOrQuit] Qwaits; 
  // // variables to record action values 
  matrix[nWaitOrQuit, N] Qwaits_ = rep_matrix(0, nWaitOrQuit, N);
  vector[N] V0_ = rep_vector(0, N);
  // // expected return 
  real G0;
  
  // initialize action values 
  //// the initial value of t = 0 
  V0 = V0_ini; 
  // the initial waiting value delines with elapsed time 
  // and the prior parameter determines at which step it falls below V0
  for(i in 1 : nWaitOrQuit){
    Qwaits[i] = - tWaits[i] * 0.1 + eta + V0;
  }
  
  // record initial action values
  Qwaits_[,1] = Qwaits;
  V0_[1] = V0;
 
  //loop over trials
  for(tIdx in 1 : (N - 1)){
    real T = Ts[tIdx]; // this trial ends on t = T
    int R = Rs[tIdx]; // payoff in this trial
    int nMadeAction = nMadeActions[tIdx]; // last decision point in this trial
    real LR; 
  
    // determine the learning rate 
    if(R > 0){
      LR = alpha;
    }else{
      LR = alphaU;
    }
    // update Qwaits towards the discounted returns
    for(i in 1 : nMadeAction){
      real t = tWaits[i]; // time for this decision points 
      real Gt = gamma ^ (T - t) * (R + V0);
      Qwaits[i] = Qwaits[i] + LR * (Gt - Qwaits[i]);
    }
    
    // update V0 towards the discounted returns 
    G0 = gamma ^ (T -(-iti)) * (R + V0);
    V0 = V0 + LR * (G0 - V0);
    
    // save action values
    Qwaits_[,tIdx+1] = Qwaits;
    V0_[tIdx+1] = V0;
  }
}
model {
  // delcare variables 
  int action; 
  vector[2] actionValues; 
  // distributions for raw parameters
  raw_alpha ~ uniform(-0.5, 0.5);
  raw_rho ~ uniform(-0.5, 0.5);
  raw_tau ~ uniform(-0.5, 0.5);
  raw_gamma ~ uniform(-0.5, 0.5);
  raw_eta ~ uniform(-0.5, 0.5);
  
  // loop over trials
  for(tIdx in 1 : N){
    real T = Ts[tIdx]; // this trial ends on t = T
    int R = Rs[tIdx]; // payoff in this trial
    int nMadeAction = nMadeActions[tIdx]; // last decision point in this trial
    
    // loop over decision points
    for(i in 1 : nMadeAction){
      // the agent wait in every decision point in rewarded trials
      // and wait except for the last decision point in non-rewarded trials
      if(R == 0 && i == nMadeAction){
        action = 2; // quit
      }else{
        action = 1; // wait
      }
      // calculate the likelihood using the soft-max function
      actionValues[1] = Qwaits_[i, tIdx] * tau;
      actionValues[2] = V0_[tIdx] * tau;
      target += categorical_logit_lpmf(action | actionValues);
    } 
  }
}
generated quantities {
 // initialize variables
  vector[2] actionValues;
  int action;
  vector[nTotalAction] log_lik = rep_vector(0, nTotalAction); // trial-wise log likelihood 
  real totalLL; // total log likelihood
  int no = 1; // action index
  
  // loop over trials
  for(tIdx in 1 : N){
    real T = Ts[tIdx]; // this trial ends on t = T
    int R = Rs[tIdx]; // payoff in this trial
    int nMadeAction = nMadeActions[tIdx]; // last decision point in this trial
    
    // loop over decision points
    for(i in 1 : nMadeAction){
      if(R == 0 && i == nMadeAction){
        action = 2; // quit
      }else{
        action = 1; // wait
      }
      // calculate the likelihood using the soft-max function
      actionValues[1] = Qwaits_[i, tIdx] * tau;
      actionValues[2] = V0_[tIdx] * tau;
      log_lik[no] =categorical_logit_lpmf(action | actionValues);
      no = no + 1;
    }
  }
  // calculate total log likelihood
  totalLL =sum(log_lik);
}
