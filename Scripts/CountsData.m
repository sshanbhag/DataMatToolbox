%------------------------------------------------------------------------
% CountsData
%------------------------------------------------------------------------
% examines spike count data
%	Files:
% 		Counts_Syllables_0-100.csv
% 		Counts_Syllables_0-200.csv
% 		Counts_Syllables_0-400.csv
% 		Counts_Syllables_0-800.csv
% 		Counts_Strings_0-100.csv
% 		Counts_Strings_0-200.csv
% 		Counts_Strings_0-400.csv
% 		Counts_Strings_0-800.csv
% 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% sharad shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created 19 April, 2013;
%
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% set path to toolboxes
%------------------------------------------------------------------------
%------------------------------------------------------------------------

%------------------------------------------------
% STA toolkit
%------------------------------------------------
tmp = which('metric');
if isempty(tmp)
	fprintf('%s: setting up STAkit paths...', mfilename)
	tmppath = '/Users/sshanbhag/Work/Code/Matlab/stable/Toolbox/STAkitMAC';
	if exist(tmppath, 'dir')
		addpath(genpath(tmppath))
	else
		error('%s: path to STAkit (SpikeTrainAnalysis Toolbox) not found', ...
						mfilename);
	end
	fprintf('... done \n\n')
end
clear tmp tmppath;
%------------------------------------------------
% Information Breakdown Toolbox
%------------------------------------------------
tmp = which('information');
if isempty(tmp)
	tmppath = '/Users/sshanbhag/Work/Code/Matlab/stable/Toolbox/InfoBreakdownToolbox';
	fprintf('%s: setting up InfoBreakdownToolbox paths...', mfilename)
	if exist(tmppath, 'dir')
		addpath(genpath(tmppath))
	else
		error('%s: path to InfoBreakdownToolbox not found', mfilename);
	end
	fprintf('... done \n\n')
end
clear tmp tmppath;

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% settings/options
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

Unit = 1;
Nbootstrap = 10;		% # of shuffled tests

%------------------------------------------------
% STA
%------------------------------------------------
q = [0 2.^(-1:5)];
q = q(1:2);
z = -2;

% settings for STAtoolkit
% possible settings: {'plugin', 'tpmc', 'jack', 'ma', 'bub', 'chaoshen', 'ww',
% 'nsb'}
opts.entropy_estimation_method = {'tpmc'};
opts.shift_cost = q;
opts.start_time = 0;
opts.end_time = .8;
opts.metric_family = 0;
opts.clustering_exponent = z;	
opts.parallel = 1;
opts.unoccupied_bins_strategy = -1;
opts.possible_words = 'recommended';


%------------------------------------------------------------------------
%------------------------------------------------------------------------
% set path to data
%------------------------------------------------------------------------
%------------------------------------------------------------------------
datapath = '/Users/sshanbhag/Work/Data/DataWave/batmat/AllSpkTrains_AllUnits';
stringfiles = {	'Counts_Strings_0-100.csv', ...
						'Counts_Strings_0-200.csv', ...
						'Counts_Strings_0-400.csv', ...
						'Counts_Strings_0-800.csv'	};
syllfiles = {	'Counts_Syllables_0-100.csv', ...
					'Counts_Syllables_0-200.csv', ...
					'Counts_Syllables_0-400.csv', ...
					'Counts_Syllables_0-800.csv'	};

%-------------------------------------------------------------------------
% load spike count data 
%-------------------------------------------------------------------------

% string data
Nstringdatafiles = length(stringfiles);
stringdata = cell(Nstringdatafiles, 1);
for n = 1:Nstringdatafiles
	stringdata{n} = csvread(fullfile(datapath, stringfiles{n}));
end
% syllable data
Nsylldatafiles = length(syllfiles);
sylldata = cell(Nsylldatafiles, 1);
for n = 1:Nsylldatafiles
	sylldata{n} = csvread(fullfile(datapath, syllfiles{n}));
end

%-------------------------------------------------------------------------
% select data
%-------------------------------------------------------------------------
data = sylldata{4};

%-------------------------------------------------------------------------
% some values from data
%-------------------------------------------------------------------------
% # of trials, stimuli, units
% depends on test type
[nRows, Trials] = size(data);
Conditions = 13; %13 for Type 1 (syllables) ; 8 for strings
nUnits = nRows/Conditions;

return
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%% compare different methods for calculating information
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------

for U = 1:1
	fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
	fprintf('\tUNIT %d\n', U);
	fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')

	%-------------------------------------------------------------------------
	% extract the raw spiketrains for the current unit
	%-------------------------------------------------------------------------
	SpkTrains = AllSpkTrains((U*Conditions+1 - Conditions):(U*Conditions),:);
	
	%-------------------------------------------------------------------------
	% first, use original method from SpikeTrain2
	%-------------------------------------------------------------------------
	fprintf('\n');
	fprintf('Original method:\n')
	clear Ddistance Dspike Dave Dsum D Hvalue
	% preallocate things
	OrigResults(U).Dconfuse = cell(Ncost, 1);
	OrigResults(U).Hvalue = zeros(Ncost, 1);	
	% loop through cost values	
	tic
	for Cindex = 1:Ncost
		Cost = q(Cindex);
		fprintf('q = %f\n', Cost);
		% allocate Distance matrix
		Ddistance = zeros(Conditions, Conditions);
		Dspike = zeros(Conditions, Trials);
		% loop through condition to use as "reference"
		for C1 = 1:Conditions
			% loop through trials
			for T1 = 1:Trials
				% for each pair of conditions (stimuli), compute spike
				% distance metric for all trials
				for C2 = 1:Conditions
					for T2=1:Trials
						Dspike(C2,T2) = (spkd(SpkTrains{C1,T1}, SpkTrains{C2,T2}, Cost))^z;
					end;
				end;

				% set distance from each train to itself to zero
				Dspike(C1,T1) = 0;   
				% find infinite values and set to 0
				[a,b] = find(Dspike==Inf);
				Dspike(a,b) = 0;

				% compute mean
				Dsum = sum(Dspike,2);
				% adjust for situations in which two conditions are identical
				for C3 = 1:Conditions,
					if C3 == C1
						Dave(C3) = (Dsum(C3)/(Trials-1))^(1/z);
					else
						Dave(C3) = (Dsum(C3)/Trials)^(1/z);
					end;
				end;
				% trap Inf values and set to 0
				[a, b] = find(Dave == Inf);
				Dave(a, b) = 0;
				% find indices of minimum avg distance
				[c ,d] = min(Dave);
				Ddistance(C1,d) = Ddistance(C1,d)+1;
				% clear things
				D(:, C1) = Dave';
				clear Dspike Dsum Dave;
			end;
		end;
		% compute some stats from Ddistance (confusion) matrix
		OrigMethod(U).Dconfuse{Cindex} = Ddistance;
		OrigMethod(U).Hvalue(Cindex) = Hstat(Ddistance);
		fprintf('\tHvalue: %f\n', OrigMethod(U).Hvalue(Cindex));
		fprintf('\tDdistance = \n');
		disp(Ddistance);
	end	% END cost loop
	orig_time = toc;
	fprintf('original method took %f seconds\n', orig_time);
	fprintf('\n');
	
	
	%-------------------------------------------------------------------------
	% bootstrap the shuffled trains
	%-------------------------------------------------------------------------
	fprintf('\n');
	fprintf('Original method, shuffled data:\n')
	clear Ddistance Dspike Dave Dsum D Hvalue
	% preallocate things
	OrigShuf(U).Dconfuse = cell(Ncost, Nbootstrap);
	OrigShuf(U).Hvalue = zeros(Ncost, Nbootstrap);	
	% loop through cost values	
	for Cindex = 1:Ncost
		Cost = q(Cindex);
		fprintf('q = %f\n', Cost);
		
		for B = 1:Nbootstrap
			% shuffle trains
			ShufTrains = shufflecell(SpkTrains);
			% allocate Distance matrix
			Ddistance = zeros(Conditions, Conditions);
			Dspike = zeros(Conditions, Trials);
			% loop through condition to use as "reference"
			for C1 = 1:Conditions
				% loop through trials
				for T1 = 1:Trials
					% for each pair of conditions (stimuli), compute spike
					% distance metric for all trials
					for C2 = 1:Conditions
						for T2=1:Trials
							Dspike(C2,T2) = (spkd(ShufTrains{C1,T1}, ShufTrains{C2,T2}, Cost))^z;
						end;
					end;

					% set distance from each train to itself to zero
					Dspike(C1,T1) = 0;   
					% find infinite values and set to 0
					[a,b] = find(Dspike==Inf);
					Dspike(a,b) = 0;

					% compute mean
					Dsum = sum(Dspike,2);
					% adjust for situations in which two conditions are identical
					for C3 = 1:Conditions,
						if C3 == C1
							Dave(C3) = (Dsum(C3)/(Trials-1))^(1/z);
						else
							Dave(C3) = (Dsum(C3)/Trials)^(1/z);
						end;
					end;
					% trap Inf values and set to 0
					[a, b] = find(Dave == Inf);
					Dave(a, b) = 0;
					% find indices of minimum avg distance
					[c ,d] = min(Dave);
					Ddistance(C1,d) = Ddistance(C1,d)+1;
					% clear things
					D(:, C1) = Dave';
					clear Dspike Dsum Dave;
				end;
			end;
			% compute some stats from Ddistance (confusion) matrix
			OrigShuf(U).Dconfuse{Cindex, B} = Ddistance;
			OrigShuf(U).Hvalue(Cindex, B) = Hstat(Ddistance);
		end	% END Bootstrap loop
		fprintf('\tHvalue (mean): %f\n', mean(OrigShuf(U).Hvalue(Cindex, :)));
		fprintf('\tHvalue (std. dev): %f\n', std(OrigShuf(U).Hvalue(Cindex, :)));
	end	% END cost loop
	fprintf('\n');

	
	
	

	%-------------------------------------------------------------------------
	% then, use STA
	%-------------------------------------------------------------------------
	fprintf('STA toolkit:\n')
	% settings for STAtoolkit
	% possible settings: {'plugin', 'tpmc', 'jack', 'ma', 'bub', 'chaoshen', 'ww',
	% 'nsb'}
	opts.entropy_estimation_method = {'tpmc'};
	opts.shift_cost = q;
	opts.start_time = 0;
	opts.end_time = .8;
	opts.metric_family = 0;
	opts.clustering_exponent = z;	
	opts.parallel = 1;
	opts.unoccupied_bins_strategy = -1;
	opts.possible_words = 'recommended';
	Ncost = length(opts.shift_cost);

	% build input using spikes2sta function
	Sin = spikes2sta(SpkTrains, 'timescale', 1, 'timeresolution', 1, ...
								'start_time', opts.start_time, 'end_time', opts.end_time);
	% call STAtoolkit function metric_shuf (same as metric, but also does shuffled
	% data procedure to check random floor
	tic
	[STAresults{U}, STAshuf{U}] = metric_shuf(Sin, opts, Nshuf);
	% pull out value for different q values
	for Cindex = 1:Ncost
		fprintf('----------------------\n')
		fprintf('\t q = %.2f\n', opts.shift_cost(Cindex));
		fprintf('----------------------\n')	
		fprintf('\tH_sta: %f\n', STAresults{U}(Cindex).table.information.value);
		fprintf('\tConfusion Matrix:\n')
		% NOTE that confusion matrix from STA toolkit is TRANSPOSE of from from long
		% (slow) method!!!
		disp(STAresults{U}(Cindex).cm');

		% get all the shuffled information bits
		hshuf = zeros(Nshuf, 1);
		for n = 1:Nshuf
			hshuf(n) = STAshuf{U}(Cindex, n).table.information.value;
		end		
		fprintf('\tShuffled Data:\n');
		fprintf('\t\tmean:\t%f\n', mean(hshuf));
		fprintf('\t\tstd:\t%f\n', std(hshuf));
		fprintf('\n')
	end
	fprintf('\n')
	sta_time = toc;
	fprintf('STAtoolkit method took %f seconds\n', sta_time);

end	% END U loop