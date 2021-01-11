# our reinfocement learning generative models simulate adapative persistence behavior as wait-or-quit choices 
# QL1: Q-learning model with a single learning rate
# QL2: Q-learning model with separate learning rates for rewards and non-rewards
# RL1: R-learning model with a single learning rate 
# RL2: R-learning model with separate learning rates for rewards and non-rewards

# inputs:
# paras: learning parameters
# condition_: HP or LP
# scheduledWait_: trial-wise delay

# outputs
# trialNum : [nTrialx1 int] 1 : nTrial
# condition : [nTrialx1 factor] from inputs
# scheduledWait : [nTrialx1 num] from inputs 
# trialEarnings : [nTrialx1 int] payment for each trial, either 10 or 0
# timeWaited : [nTrialx1 num] how long the agent waits after the iti in each trial 
# sellTime : [nTrialx1 num]  when the agent sells the token on each trial 
# Qwaits_ : [20/40 x nTrial num] value of waiting at each second in each trial
# V_ : [nTrialx1 num] value of entering a pre-trial iti, namely t = 0

RL1 = function(paras, condition, duration, normResults, delays = c()){
  source("subFxs/taskFxs.R")
  
  # check whether delays have been given 
  exist = length(delays) > 0
  
  # default settings 
  iti = 2
  
  # normative analysis 
  optimRewardRates = normResults$optimRewardRates
  optimWaitThresholds = normResults$optimWaitThresholds
  
  # learning parameters
  alpha = paras[1]; tau = paras[2];  prior = paras[3]; beta = paras[4];
  
  # num of trials
  nTrialMax = ceiling(duration / iti)
  # duration of a sampling interval 
  stepSec = 1
  # max delay duration 
  delayMax = ifelse(condition == "HP", delayMaxs[1], delayMaxs[2])
    
  # initialize action values 
  V = 0 # state value for t = 0
  rewardRate = mean(unlist(optimRewardRates))
  tWaits = seq(0, delayMax, by = stepSec) # decision points 
  tMax = max(tWaits) #  time point for the last decision point
  Qwaits = -0.1 * (tWaits) + prior + V # value of waiting at each decision points
  
  # initialize output variables
  Gs_ = matrix(0, length(tWaits), nTrialMax)
  Qwaits_ = matrix(NA, length(tWaits), nTrialMax); Qwaits_[,1] = Qwaits 
  V_ = vector(length = nTrialMax); V_[1] = V
  scheduledWait_ =  rep(NA, nTrialMax)
  trialEarnings_ = rep(NA, nTrialMax)
  timeWaited_ = rep(NA, nTrialMax)
  sellTime_ = rep(NA, nTrialMax)
  
  # track elpased time from the beginning of the task 
  elapsedTime = -iti # the first trial doesn't have a pre-trial ITI 
  
  # loop over trials
  tIdx = 1
  while(elapsedTime < duration){
    # current scheduledWait 
    if(exist){
      scheduledWait = delays[tIdx]
    }else{
      scheduledWait = drawSample(condition)
    }
    scheduledWait_[tIdx] = scheduledWait
      
    # sample at a temporal resolution of 1 sec until a trial ends
    t = -2
    while(t <= tMax){
      if(t >= 0){
        # decide whether to wait or quit
        pWait = 1 / sum(1  + exp((V - Qwaits[tWaits == t])* tau))
        action = ifelse(runif(1) < pWait, 'wait', 'quit')
        
        # if a reward occurs and the agent is still waiting, the agent gets the reward
        tokenMature = (scheduledWait >= t) & (scheduledWait < (t + stepSec)) # whether the token matures before the next decision point
        getToken = (action == 'wait' && tokenMature) # whether the agent obtains the matured token
        
        # a trial ends if the agent obtains the matured token or quits. 
        # if the trial ends,return to t = 0. Otherwise, proceed to t + 1.
        isTerminal = (getToken || action == "quit")
        if(isTerminal){
          # update trial-wise variables 
          T =  ifelse(getToken, scheduledWait, t) # if quit, time ends at t
          timeWaited =  T # how long the agent waits since the token appears
          trialEarnings = ifelse(getToken, tokenValue, 0) 
          elapsedTime = elapsedTime + timeWaited + iti
          sellTime = elapsedTime # elapsed task time when the agent sells the token
          # record trial-wise variables
          trialEarnings_[tIdx] = trialEarnings
          timeWaited_[tIdx] = timeWaited
          sellTime_[tIdx] = sellTime
          break
        }
      }
      t = t + stepSec
    }
    
    # when the trial endes, update value functions for all time points before T in the trial

    # calculate expected returns for t >= 2
    Gts = trialEarnings - rewardRate * (T - tWaits) + V
    # only update value functions before time t = T
    updateFilter = tWaits <= T 
    Gs_[updateFilter,tIdx] = Gts[updateFilter]
    # update Qwaits
    Qwaits[updateFilter] = Qwaits[updateFilter] + alpha * (Gts[updateFilter] - Qwaits[updateFilter])
    
    # calculate expected returns for t == 0 and update V
    Gt =  trialEarnings - rewardRate * (T - (-2)) + V
    delta = Gt - V
    V = V + alpha * delta
    rewardRate = rewardRate + beta * delta
    
    # record updated values
    Qwaits_[,tIdx + 1] = Qwaits
    V_[tIdx + 1] = V
    
    # proceed to the next trial
    tIdx = tIdx + 1
  } # end of the loop over trials
  
  # return outputs
  nTrial = min(which(is.na(sellTime_))) - 1
  outputs = list( 
    "trialNum" = 1 : nTrial, 
    "condition" = rep(condition, nTrial),
    "trialEarnings" = trialEarnings_[1 : nTrial], 
    "timeWaited" = timeWaited_[1 : nTrial],
    "sellTime" = sellTime_[1 : nTrial],
    "scheduledWait" = scheduledWait_[1 : nTrial],
    "Qwaits_" = Qwaits_[, 1 : nTrial], 
    "Gs_" = Gs_[, 1 : nTrial], 
    "V_" = V_[1 : nTrial]
  )
  return(outputs)
}