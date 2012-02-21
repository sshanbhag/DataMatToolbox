%------------------------------------------------------------------------
%------------------------------------------------------------------------
% load data
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% load(fullfile(inpath, 'statsdata.mat'));
load('statsdata.mat')

%------------------------------------------------------------------------
% analysis parameters
%------------------------------------------------------------------------
% power transformation for averaging, usually -2
% not used in original SpikeTrain.m function - added back here 
% using jmatrix.m as model
z = -2;


%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Analyze LFH Data
%------------------------------------------------------------------------
%------------------------------------------------------------------------

Nunits = length(Mlfhstruct{1});
Nconditions = length(Mlfhstruct);
ConditionList = 1:Nconditions

uIndex = 1;
% assign data to M struct for simplicity
for condition = 1:Nconditions
	M(condition) = Mlfhstruct{condition}(uIndex);
end

fprintf('File %s, uIndex: %d\n', M(1).file, uIndex);
fprintf('  Probe %d, Cluster %d\n', M(1).probe, M(1).cluster);

Ntrials = M(1).ntrials;
Ntrials = 4;


CostVals = [0, 0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1];

CostVals = [0.6];
opts.shift_cost = 0.001 .* CostVals;
opts.metric_family=0;

Times = cell(Ntrials * Nconditions, 1);
Labels = cell(Ntrials * Nconditions, 1);
Categories = zeros(Ntrials * Nconditions, 1);

tindex = 0;
for c = 1:Nconditions
	for n = 1:Ntrials
		tindex = tindex + 1;
		Times{tindex} = 0.001 .* M(c).spikes{n};
		Labels{tindex} = ones(size(M(c).spikes{n}));
		Categories(tindex) = c;
	end
end

D = metricdist(1, Times, Labels, opts);

CM = metricclust(D, int32(Categories), 3);

%{
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
%}



%------------------------------------------------------------------------
% loop units
%------------------------------------------------------------------------
Ncost = length(CostVals);

uIndex = 1;
	
% assign data to M struct for simplicity
for condition = 1:Nconditions
	M(condition) = Mlfhstruct{condition}(uIndex);
end

fprintf('File %s, uIndex: %d\n', M(1).file, uIndex);
fprintf('  Probe %d, Cluster %d\n', M(1).probe, M(1).cluster);

% loop through cost values
for Cindex = 1:Ncost
	% get current cost value from list
	Cost = CostVals(Cindex);
	fprintf('\tCost = %f\n', Cost);

	% initialize Ddistance array
	Ddistance = zeros(Nconditions, Nconditions);

	% initialize Dspike
	Dspike = zeros(Nconditions*Ntrials);
	dind = cell(size(Dspike));

	% initialize rowIndex (for indexing into rows of Dspike)
	rowIndex = 0;
	colIndex = 0;
	
	% C1 loop for the number of conditions
	for C1 = 1:Nconditions
		
		% T1 loop for the number of trials;
		for T1 = 1:Ntrials
		
			% increment row index
			rowIndex = rowIndex + 1;
	
			% for each pair of conditions (stimuli), compute spike
			% distance metric for all trials (except for this
			% trial/condition pair with itself)
			for C2 = 1:Nconditions
				for T2=1:Ntrials
					% increment colIndex
					colIndex = colIndex + 1;
					tmp = spkd(M(C1).spikes{T1}, M(C2).spikes{T2}, Cost);
					Dspike(rowIndex,colIndex) = tmp^z;
					dind{rowIndex, colIndex} = [C1, T1; C2, T2];
				end;	% END T2 LOOP
			end;	% END C2 LOOP
			% reset column index ( for indexing columns of Dspike)
			colIndex = 0;
		end	%	END C1 LOOP
	end	%	END T1 LOOP

	% find values of Dspike that are 'Inf' and set them to zero
 	[zeroind] = find(Dspike==Inf);
 	Dspike(zeroind) = 0;

	%{
	% compute overall sum across trials (rows of Dspike matrix)
	Dsum = sum(Dspike,2);

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
	Ddistance(C1, d) = Ddistance(C1, d) + 1;

	% compute some stats from Ddistance (confusion) matrix
	[Hvalue, Hrand, Hprob] = Hstat(Ddistance);

	% ExsName= strcat(Filename.name,'-Cell-', num2str(Units))
	% xlswrite(ExsName, Ddistance)
	fprintf('\tHvalue: %f\n', Hvalue);
	fprintf('\tHrand: %f\n', Hrand);
	fprintf('\tHprob: %f\n', Hprob);
	fprintf('\tDdistance = \n');
	disp(Ddistance);

	% store data
	Dlfh(uIndex).Ddistance{Cindex} = Ddistance;
	Dlfh(uIndex).Hvalue(Cindex) = Hvalue;
	Dlfh(uIndex).Hrand(Cindex) = Hrand;
	Dlfh(uIndex).Hprob(Cindex) = Hprob;

	clear Ddistance SpkTrains;
	
	%}
end


