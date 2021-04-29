
nExp = 3;

% compare the four variants of the RL model
EF4_ = zeros([3, 4]); % expeted model frequencies 
PXP4_ = zeros([3, 4]); % protected prob of exceedance 
for i = 1 : nExp
	fileDir = sprintf('../../genData/wtw_exp%d/waic.csv', i);
    % fileDir = sprintf('waic_exp%d.csv', i);
	waic = csvread(fileDir);
	
	% transfer into the loss domain 
	L = -waic';

	% exclude participants with disconvergent model fitting results 
	L = L(1:4, all(L(7:10,:)));

	% run the group-level Bayesian model selection, Stephan et al., 2009 
	[posterior,out] = VBA_groupBMC(L);
	f  = out.Ef;
	EP = out.ep;
	EF4_(i, :)  = out.Ef;
	PXP4_(i, :) = out.ep;
end 
csvwrite('../../genData/EF4.csv', EF4_, 3, 4);
csvwrite('../../genData/PXP4.csv', PXP4_, 3, 4);

% add the two non-learning benchmark models 
EF6_ = zeros([3, 6]); % expeted model frequencies 
PXP6_ = zeros([3, 6]); % protected prob of exceedance 
for i = 1 : nExp
	fileDir = sprintf('../../genData/wtw_exp%d/waic.csv', i);
	waic = csvread(fileDir);
	
	% transfer into the loss domain 
	L = -waic';

	% exclude participants with disconvergent model fitting results 
	L = L(1:6, all(L(7:12,:)));

	% run the group-level Bayesian model selection, Stephan et al., 2009 
	[posterior,out] = VBA_groupBMC(L);
	f  = out.Ef;
	EP = out.ep;
	EF6_(i, :)  = out.Ef;
	PXP6_(i, :) = out.ep;
end 
csvwrite('../../genData/EF6.csv', EF6_, 3, 6);
csvwrite('../../genData/PXP6.csv', PXP6_, 3, 6);