%------------------------------------------------------------------------
% mextest
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Used to evaluate STAtoolkit from Victor et al.
%------------------------------------------------------------------------
% 3 April, 2013, SJS
%------------------------------------------------------------------------
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% clear and cleanup
%------------------------------------------------------------------------
%------------------------------------------------------------------------
clear all

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% set path to toolbox
%------------------------------------------------------------------------
%------------------------------------------------------------------------
tmp = which('metric');
if isempty(tmp)
	fprintf('%s: setting up Spike paths...', mfilename)
	addpath(genpath('/Users/sshanbhag/Work/Code/Matlab/dev/Toolboxes/STAkitMAC'))
	fprintf('... done \n\n')
end
clear tmp;

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% load datums
%------------------------------------------------------------------------
%------------------------------------------------------------------------
datapath = '/Users/sshanbhag/Work/Data/DataWave/batmat/AllSpkTrains_AllUnits'
datafile = 'AllSpkTrains_AllUnits0to200_strings.mat';
load(fullfile(datapath, datafile));

% # of trials, stimuli, units depends on test type
% strings = 8 stimuli
% sylls = 13 stimuli
[nTrials, nReps]=size(AllSpkTrains);
nStimuli=8; %13 for Type 1 (syllables) ; 8 for strings
nUnits=nTrials/nStimuli;

Nbootstrap = 5;

return

for unit = 1:10

	SpkTrains = AllSpkTrains((unit*nStimuli+1 - nStimuli):(unit*nStimuli),:);


	spikes = spikes2sta(SpkTrains, 'timescale', 1, 'timeresolution', 1)

	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% The options and parameters for metric:
	%------------------------------------------------------------------------
	% 	OPTS.start_time:		The start time of the analysis window. 
	% 									The default is the maximum of all 
	% 									of the start times in X. 
	% 
	% 	OPTS.end_time:			The end time of the analysis window. 
	% 									The default is the minimum of all of the 
	% 									end times in X.
	% 
	% 	OPTS.shift_cost:		The cost of shifting a spike per unit time
	% 								relative to inserting or deleting a spike. 
	% 								This option may be a vector of such values. 
	% 									The default is  1/(end_time-start_time).
	% 
	% 	OPTS.label_cost:		This applies only to data sets with simultaneously 
	% 								recorded spike trains. It is the cost of
	% 								altering a spike's label, and may range from 0 to 2. 
	% 								This option may be a vector of such values whose length 
	% 								is equal to OPTS.shift_cost. 
	% 									The default is 0.
	% 
	% 	OPTS.metric_family:	Selects the metric to be used.
	% 								  The default value is 0.
	% 	  OPTS.metric_family=0:		Uses D^spike metric.
	% 	  OPTS.metric_family=1:		Uses D^interval metric. This is
	% 										only applicable to single-site data.
	% 
	% 	OPTS.parallel:			Selects which algorithm version to use.
	% 									  The default value is 0 if OPTS.shift_cost 
	% 									  has one element	and 1 if OPTS.shift_cost 
	% 									  has multiple elements.
	% 	  OPTS.parallel=0:	Computes distances for a single shift_cost,
	% 								label_cost pairs at a time.
	% 	  OPTS.parallel=1:	Uses an algorithm that computes the
	% 								distances for all shift_cost,label_cost pairs
	% 								concurrently. When many parameters sets are being
	% 								analyzed, this method can provide considerable
	% 								computational savings. 
	% 
	% 	OPTS.clustering_exponent:	A constant that controls the
	% 									  clustering. Negative values emphasize 
	% 									  smaller distances and positive values 
	% 									  emphasize larger distances. 
	% 										The default is -2.
	%------------------------------------------------------------------------
	opts.entropy_estimation_method = {'plugin'};
	opts.shift_cost = [0 2.^(-2:2)];
	opts.start_time = 0;
	opts.end_time = 1;
	opts.metric_family = 0;
	opts.clustering_exponent = 2;
	opts.parallel = 1;
	opts.unoccupied_bins_strategy = -1;
	Ncost = length(opts.shift_cost);

	% a = staread('/Users/sshanbhag/Work/Code/Matlab/dev/Toolboxes/STAkitMAC/data/taste.stam')
	% [aout, optsout] = metric(a, opts)

	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% use toolkit to obtain info measures
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% spike distance metric
	% [xout, optsout] = metric(x, opts);
	% shuffled metric
	[x, s, soptsout] = metric_shuf(spikes, opts, Nbootstrap);

	for c = 1:Ncost
	  xbits(c) = x(c).table.information.value;
	end

	sbits_all = zeros(Ncost, Nbootstrap);
	for c = 1:Ncost
		for b = 1:Nbootstrap
			sbits_all(c, b) = s(c, b).table.information.value;
		end
	end
	sbits_mean = mean(sbits_all, 2);
	sbits_sd = std(sbits_all, 0, 2);

	[max_info, max_indx] = max(xbits);

	figure
	subplot(211)
	plot(opts.shift_cost, xbits, '.:');
	xlabel('Temporal precision (1/sec)');
	ylabel('Information (bits)');
	hold on
		errorbar(opts.shift_cost, sbits_mean, sbits_sd, 'r.:')
	hold off
	title({sprintf('UNIT %d', unit), 'STA method'})
	drawnow
	
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% long way to do it
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------

	%------------------------------------------------------------------------
	%% loop through cost values
	%------------------------------------------------------------------------
	for Cindex = 1:Ncost

		% get current cost value from list
		Cost = opts.shift_cost(Cindex);
		fprintf('Cost = %f\n', Cost);

		% initialize Ddistance array
		Ddistance= zeros(nStimuli, nStimuli);

		%------------------------------------------------------------------------
		% compute UNSHUFFLED metric
		%------------------------------------------------------------------------
		Ddistance = confusematrix(SpkTrains, Cost, opts.clustering_exponent);		
		% compute some stats from Ddistance (confusion) matrix
		Hvalue = Hstat(Ddistance);

		fprintf('Unit %d\n', unit)
		fprintf('\tHvalue: %f\n', Hvalue);
		disp(Ddistance);
		Ddata(unit).Ddistance{Cindex} = Ddistance;
		Ddata(unit).Hvalue(Cindex) = Hvalue;

		%------------------------------------------------------------------------
		% compute SHUFFLED metric
		%------------------------------------------------------------------------
		fprintf('Computing bootstrap/shuffled values...');
		Hb = zeros(Nbootstrap, 1);
		for B = 1:Nbootstrap
			% shuffle spike trains across Stimuli and reps
			ShuffleTrains = shufflecell(SpkTrains);
			Dshuf = confusematrix(ShuffleTrains, Cost, opts.clustering_exponent);
			% compute some stats from Ddistance (confusion) matrix
			Hb(B) = Hstat(Dshuf);
		end	% END B loop
		Ddata(unit).Hb{Cindex} = Hb;
		fprintf('... done!\n')
		fprintf('\t\tmean(Hb) = %f\n', mean(Hb));

	end	% END Cindex loop (cost values)

	%------------------------------------------------------------------------
	%% plots
	%------------------------------------------------------------------------
	subplot(2, 1, 2)
	plot(opts.shift_cost, Ddata(unit).Hvalue, 'b.-')
	for c = 1:Ncost
		Hbmean(c) = mean(Ddata(unit).Hb{c});
		Hbstd(c) = std(Ddata(unit).Hb{c});
	end
	hold on
		errorbar(opts.shift_cost, Hbmean, Hbstd, 'r.-')
	hold off
	xlabel('Temporal precision (1/sec)');
	ylabel('Information (bits)');
	title('long method')

end
	