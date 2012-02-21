function [Dclust, H] = metric_cluster(Dspike, Costs, Cats, Z)
%-----------------------------------------------------------------------------
% [Dclust, H] = metric_cluster(Dspike, Costs, Cats, Z)
%-----------------------------------------------------------------------------
% 
% Description
% 
%-----------------------------------------------------------------------------
% Input Arguments:
% 	Input		input info
% 
% Output Arguments:
% 	Output	output info
%
%-----------------------------------------------------------------------------
% See also: 
%-----------------------------------------------------------------------------

%-----------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%-----------------------------------------------------------------------------
% Created: 21 February, 2012 (SJS)
%
% Revisions:
%-----------------------------------------------------------------------------
% TO DO:
%-----------------------------------------------------------------------------


%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Check inputs
%------------------------------------------------------------------------
%------------------------------------------------------------------------

if ~exist('Z', 'var')
	Z = -2;
	fprintf('%s: using default Z value (%d)\n', mfilename, Z);
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% extract some needed information and initialize things
%------------------------------------------------------------------------
%------------------------------------------------------------------------
Ncosts = length(Costs);

ConditionList = unique(Cats);
Nconditions = length(ConditionList);
% initialize list of indices into D matrix by condition
conditionIndices = cell(Nconditions, 1);
Ntrials = zeros(Nconditions, 1);
for c = 1:Nconditions
	conditionIndices{c} = (Cats == c);
	Ntrials(c) = nnz(conditionIndices{c});
end
TotalTrials = sum(Ntrials);


% initialize Dclust array
Dclust = zeros(Nconditions, Nconditions, Ncosts);


%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Cluster analysis routine
%------------------------------------------------------------------------
%------------------------------------------------------------------------

% apply exponent to data
Dspike_exp = Dspike.^Z; 

% find values of Dspike_exp that are 'Inf' and set them to zero
Dspike_exp(Dspike_exp==Inf) = 0;

H = repmat(struct('value', [], 'rand', [], 'prob', []), Ncosts, 1);

% loop through cost values
for Cindex = 1:Ncosts
	% get current cost value from list
	fprintf('\tCost = %f\n', Costs(Cindex));
	
	% get adjusted distances for current cost value
	tmpD = Dspike_exp(:, :, Cindex);
	tmpD_bycondition = cell(Nconditions, 1);

	% loop for the TOTAL number of trials (# rows in Dspike)
	for T = 1:TotalTrials
		
		% get the current trials' condition from the Category list
		testCondition = Cats(T);
		% get the other conditions from the Conditions
		otherConditions = ConditionList(testCondition ~= ConditionList);
		
		%-----------------------------------------------------------------------
		%-----------------------------------------------------------------------
		% get average distance from this trial to other trials
		%-----------------------------------------------------------------------
		%-----------------------------------------------------------------------
		
		%-----------------------------------------------------------------------
		% first, break up data by condition	
		%-----------------------------------------------------------------------
		% current condition
		tmpD_bycondition{testCondition} = tmpD(T, conditionIndices{testCondition});
		% other conditions
		for c = 1:length(otherConditions)
			tmpind = otherConditions(c);
			tmpD_bycondition{tmpind} = tmpD(T, conditionIndices{tmpind});
		end
		
		%-----------------------------------------------------------------------
		% now, compute sum across trials
		%-----------------------------------------------------------------------
		tmpSum = zeros(1, Nconditions);
		for c = 1:Nconditions
			tmpSum(c) = sum(tmpD_bycondition{c});
		end
		
		%-----------------------------------------------------------------------
		% compute averages (dropping current trial from computation)
		%-----------------------------------------------------------------------
		tmpAve = zeros(1, Nconditions);
		tmpAve(testCondition) = tmpSum(testCondition) ./ ...
												(Ntrials(testCondition) - 1);
		for c = 1:length(otherConditions)
			tmpind = otherConditions(c);
			tmpAve(tmpind) = tmpSum(tmpind) ./ Ntrials(tmpind);
		end

		%-----------------------------------------------------------------------
		% undo the exponent adjustment
		%-----------------------------------------------------------------------
		tmpAve = tmpAve.^(1/Z);
		
		% find value and indices of minimum avg distance
		[~, mind] = min(tmpAve);
		% increment count for this overall condition
		Dclust(testCondition, mind, Cindex) = Dclust(testCondition, mind, Cindex) + 1;
	end

	% compute some stats from Ddistance (confusion) matrix
	[H(Cindex).value, H(Cindex).rand, H(Cindex).prob] = Hstat(Dclust(:, :, Cindex));

	% ExsName= strcat(Filename.name,'-Cell-', num2str(Units))
	% xlswrite(ExsName, Ddistance)
	fprintf('\tHvalue: %f\n', H(Cindex).value);
	fprintf('\tHrand: %f\n', H(Cindex).rand);
	fprintf('\tHprob: %f\n', H(Cindex).prob);
	fprintf('\tDclust = \n');
	disp(Dclust(:, :, Cindex));
end
