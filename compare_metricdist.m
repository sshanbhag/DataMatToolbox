%------------------------------------------------------------------------
%------------------------------------------------------------------------
% create test data
%------------------------------------------------------------------------
%------------------------------------------------------------------------
Ntrials = [4 4 4 8];
TotalTrials = sum(Ntrials);
Nconditions = length(Ntrials);
ConditionList = 1:Nconditions;


Times = cell(TotalTrials, 1);
Labels = cell(TotalTrials, 1);
Categories = zeros(TotalTrials, 1);

base_spikes = 100:50:1000;
mexIndex = 0;
for c = 1:Nconditions
	for t = 1:Ntrials(c)
		Spikes{c, t} = base_spikes + mexIndex;
		mexIndex = mexIndex + 1;
		Times{mexIndex} = Spikes{c, t};
		Labels{mexIndex} = ones(size(Spikes{c, t}));
		Categories(mexIndex) = c; 
	end
end
mexIndex
TotalTrials

%------------------------------------------------------------------------
% analysis parameters
%------------------------------------------------------------------------
% power transformation for averaging, usually -2
% not used in original SpikeTrain.m function - added back here 
% using jmatrix.m as model
z = -2;

% cost values to test
CostVals = [0, 0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1];
CostVals = [.1 .2 .4 .8];
CostVals = .1;
Ncost = length(CostVals);


%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Analyze LFH Data
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% some options for metricdist() mex function
distopts.shift_cost = CostVals;
distopts.metric_family = 0;
distopts.parallel = 0;

tic
for Cindex = 1:Ncost
	distopts.shift_cost = CostVals(Cindex);
	D(:, :, Cindex) = metricdist(1, Times, Labels, distopts);
end
time_metricdist = toc


%------------------------------------------------------------------------
%------------------------------------------------------------------------
% metricclust() mexfunction that crashes and takes down Matlab... :(
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% clustopts.clustering_exponent = z;
% CM = metricclust(D, int32(Categories), int32(Nconditions), clustopts);
%------------------------------------------------------------------------
%------------------------------------------------------------------------

%{
%------------------------------------------------------------------------
METRICDIST Compute distances between sets of spike train pairs.
  D = METRICDIST(N,TIMES,LABELS,OPTS) uses the metric
  space method to compute the distances between all possible spike
  train pairs in the data set. D is a 3-D matrix where the third dimension
  corresponds to the elements of OPTS.shift_cost.

  N is the number of neurons recorded simultaneously.

  TIMES is an Px1 cell array of spike times where P is the number of
  trials in the data set. Each member of the cell array consists of
  spike times interleaved from simultaneous trials.

	each element in Times will have N X nspikes values

  LABELS is a Px1 cell array with labels of the spikes in TIMES.

		each element will have N X nspikes values identifying the units


  The options and parameters for this function are:
     OPTS.shift_cost: The cost of shifting a spike per unit time
        relative to inserting or deleting a spike. This option may
        be a vector of such values. The default is
        1/(end_time-start_time).
     OPTS.label_cost: This applies only to data sets with
        simultaneously recorded spike trains. It is the cost of
        altering a spike's label, and may range from 0 to 2. This
        option may be a vector of such values whose length is equal
        to OPTS.shift_cost. The default is 0.
     OPTS.metric_family: Selects the metric to be used.
        OPTS.metric_family=0: Uses D^spike metric.
        OPTS.metric_family=1: Uses D^interval metric. This is
           only applicable to single-site data.
        The default value is 0.
     OPTS.parallel: Selects which algorithm version to
        use.
        OPTS.parallel=0: Computes distances for a single shift_cost,
            label_cost pairs at a time.
        OPTS.parallel=1: Uses an algorithm that computes the
            distances for all shift_cost,label_cost pairs
            concurrently. When many parameters sets are being
            analyzed, this method can provide considerable
            computational savings. 
        The default value is 0 if OPTS.shift_cost has one element
            and 1 if OPTS.shift_cost has multiple elements.

  D = METRICDIST(N,TIMES,LABELS) uses the default options and
  parameters.

  [D,OPTS_USED] = METRICDIST(N,TIMES,LABELS) or [D,OPTS_USED] =
  METRICDIST(N,TIMES,LABELS,OPTS) additionally return the
  options used.
%------------------------------------------------------------------------

%------------------------------------------------------------------------
METRICCLUST Cluster spike trains based on distance matrix.
  CM = METRICCLUST(D,CATEGORIES,M,OPTS) uses a simple clustering
  method to classify spike trains based on the distance matrix
  D. CATEGORIES is a vector that gives the category indices of the
  spike trains. M is the number of categories. CM is a square matrix
  where the columns correspond to the actual classes and the rows
  correspond to the assigned classes.

  *On mac I64 systems (at least), the CATEGORIES input must be
   in int32 format.  This can be ensured by calling METRICCLUST
   as follows:  CM = metricclust(D, int32(CATEGORIES), M)

  The options and parameters for this function are:
     OPTS.clustering_exponent: A constant that controls the
        clustering. Negative values emphasize smaller distances
        and positive values emphasize larger distances. The
        default is -2.

  CM = METRICCLUST(D,CATEGORIES,M) uses the default options and parameters.

  [CM,OPTS_USED] = METRICCLUST(D,CATEGORIES,M,OPTS) or
  [CM,OPTS_USED] = METRICCLUST(D,CATEGORIES,M,OPTS) additionally
%------------------------------------------------------------------------
  return the options used.
%}


%------------------------------------------------------------------------
% loop units
%------------------------------------------------------------------------
% initialize Ddistance array
Ddistance = zeros(Nconditions, Nconditions, Ncost);
% initialize Dspike
Dspike = zeros(TotalTrials, TotalTrials, Ncost);
dind = cell(size(Dspike));

tic
% get current cost value from list
for Cindex = 1:Ncost
	Cost = CostVals(Cindex);

	% initialize rowIndex (for indexing into rows of Dspike)
	rowIndex = 0;
	colIndex = 0;

	% C1 loop for the number of conditions
	for C1 = 1:Nconditions
		% T1 loop for the number of trials;
		for T1 = 1:Ntrials(C1)

			% increment row index
			rowIndex = rowIndex + 1;

			% for each pair of conditions (stimuli), compute spike
			% distance metric for all trials (except for this
			% trial/condition pair with itself)
			for C2 = 1:Nconditions
				for T2=1:Ntrials(C2)
					% increment colIndex
					colIndex = colIndex + 1;
					Dspike(rowIndex, colIndex, Cindex) = spkd(Spikes{C1, T1}, Spikes{C2, T2}, Cost);
					dind{rowIndex, colIndex, Cindex} = [C1, T1; C2, T2];
				end;	% END T2 LOOP
			end;	% END C2 LOOP
			
			% reset column index ( for indexing columns of Dspike)
			colIndex = 0;
		end	%	END C1 LOOP
	end	%	END T1 LOOP
end
spkd_time = toc

% check difference in methods
sum(sum(sum(Dspike - D)))

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Cluster analysis routine
%------------------------------------------------------------------------
%------------------------------------------------------------------------

% apply exponent to data
Dspike_exp = Dspike.^z; 

% find values of Dspike_exp that are 'Inf' and set them to zero
zeroind = find(Dspike_exp==Inf);
Dspike_exp(zeroind) = 0;

% initialize list of indices into D matrix by condition
for c = 1:Nconditions
	conditionIndices{c} = (Categories == c);			
end

% initialize Dclust array
Dclust = zeros(Nconditions, Nconditions, Ncost);

% loop through cost values
for Cindex = 1:Ncost
	% get current cost value from list
	Cost = CostVals(Cindex);
	fprintf('\tCost = %f\n', Cost);
	
	% get adjusted distances for current cost value
	tmpD = Dspike_exp(:, :, Cindex);


	% loop for the TOTAL number of trials (# rows in Dspike)
	for T = 1:TotalTrials
		
		% get the current trials' condition from the Category list
		testCondition = Categories(T);
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
		for c = otherConditions
			tmpD_bycondition{c} = tmpD(T, conditionIndices{c});
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
		tmpAve(otherConditions) = tmpSum(otherConditions) ./ ...
												Ntrials(otherConditions);
		tmpAve
		%-----------------------------------------------------------------------
		% undo the exponent adjustment
		%-----------------------------------------------------------------------
		tmpAve = tmpAve.^(1/z);
		
		% find value and indices of minimum avg distance
		[minv, mind] = min(tmpAve);
		% increment count for this overall condition
		Dclust(testCondition, mind, Cindex) = Dclust(testCondition, mind, Cindex) + 1;
	end

	% compute some stats from Ddistance (confusion) matrix
	[Hvalue, Hrand, Hprob] = Hstat(Dclust(:, :, Cindex));

	% ExsName= strcat(Filename.name,'-Cell-', num2str(Units))
	% xlswrite(ExsName, Ddistance)
	fprintf('\tHvalue: %f\n', Hvalue);
	fprintf('\tHrand: %f\n', Hrand);
	fprintf('\tHprob: %f\n', Hprob);
	fprintf('\tDclust = \n');
	disp(Dclust);

	% store data
	V.Dclust{Cindex} = Dclust(:, :, Cindex);
	V.Hvalue(Cindex) = Hvalue;
	V.Hrand(Cindex) = Hrand;
	V.Hprob(Cindex) = Hprob;

end





% not sure if this cluster routine is viable
%{
for Cindex = 1:Ncost
	% apply exponent and compute overall sum across trials (rows of Dspike matrix)
	Dsum = sum(Dspike_exp(:, :, Cindex), 2);

	% loop through conditions, using appropriate # of trials to 
	% compute average for situation in which loop condition == current 
	% overall condition for this cost value
	Dave = zeros(Nconditions, 1);
	currentCondition = ConditionList(C1 == ConditionList);
	otherNconditions = ConditionList(C1 ~= ConditionList);
	Dave(currentCondition) = (Dsum(currentCondition)./(Ntrials-1)).^(1/z);
	Dave(otherNconditions) = (Dsum(otherNconditions)./Ntrials).^(1/z);

	% find indices of minimum avg distance
	[c ,d] = min(Dave);
	% increment count for this overall condition
	Ddistance(C1, d, Cindex) = Ddistance(C1, d, Cindex) + 1;

	% compute some stats from Ddistance (confusion) matrix
	[Hvalue, Hrand, Hprob] = Hstat(Ddistance(:, :, Cindex));

	% ExsName= strcat(Filename.name,'-Cell-', num2str(Units))
	% xlswrite(ExsName, Ddistance)
	fprintf('\tHvalue: %f\n', Hvalue);
	fprintf('\tHrand: %f\n', Hrand);
	fprintf('\tHprob: %f\n', Hprob);
	fprintf('\tDdistance = \n');
	disp(Ddistance);

	% store data
	Cdata.Hvalue(Cindex) = Hvalue;
	Cdata.Hrand(Cindex) = Hrand;
	Cdata.Hprob(Cindex) = Hprob;

end
%}
