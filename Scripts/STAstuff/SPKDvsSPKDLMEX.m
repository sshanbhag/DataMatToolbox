%------------------------------------------------------------------------
% CompareMethods_RealData
%------------------------------------------------------------------------
% compares matlab script based spike distance metric calculation with
% STAtoolkit implementation.  Test spiketrain data are from 
% stringfile = 'AllSpkTrains_AllUnits0to800_strings.mat';
% syllfile = 'AllSpkTrains_AllUnits0to800_syllables.mat';
% 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% sharad shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created 13 April, 2013;
%
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% set path to toolboxes
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% STAkitMAC
tmp = which('metric');
if isempty(tmp)
	fprintf('%s: setting up STAtoolkit paths...', mfilename)
	addpath(genpath('/Users/sshanbhag/Work/Code/Matlab/stable/Toolbox/STAkitMAC'))
	fprintf('... done \n\n')
end
clear tmp;
% SDM tools (with spkdl MEX file)
tmp = which('spkdl');
if isempty(tmp)
	fprintf('%s: setting up SDM paths...', mfilename)
	addpath(genpath('/Users/sshanbhag/Work/Code/Matlab/stable/Toolbox/SDM'))
	fprintf('... done \n\n')
end
clear tmp;

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% misc things
%------------------------------------------------------------------------
%------------------------------------------------------------------------
sepstr = '---------------------------------------------------------------';

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% set path to data
%------------------------------------------------------------------------
%------------------------------------------------------------------------
datapath = '/Users/sshanbhag/Work/Data/DataWave/batmat/AllSpkTrains_AllUnits';
stringfile = 'AllSpkTrains_AllUnits0to800_strings.mat';
syllfile = 'AllSpkTrains_AllUnits0to800_syllables.mat';

%-------------------------------------------------------------------------
% SDM settings
%-------------------------------------------------------------------------
% cost (q) in 1/seconds
q = [0 2.^(-1:5)];
q = q(1:2);
% exponent
z = -2;
Nshuf = 2;		% # of shuffled/bootstrap tests

%-------------------------------------------------------------------------
% load spike trains 
%-------------------------------------------------------------------------
% MG settings:
%load('E:\Bat - SingleCh Restrained\batmat\Converted Dobj files\ALL Dobj files2\AllSpkTrains_testnewinfo_sec_0to100_strings.mat');
%load('E:\Bat - SingleCh Restrained\batmat\Converted Dobj files\ALL Dobj files2\AllSpkTrains_TestArray_2.mat');

% SJS 
load(fullfile(datapath, syllfile));

%-------------------------------------------------------------------------
% some values from data
%-------------------------------------------------------------------------
% # of trials, stimuli, units
% depends on test type
[nRows, Trials] = size(AllSpkTrains);
Conditions = 13; %13 for Type 1 (syllables) ; 8 for strings
nUnits = nRows/Conditions;
Ncost = length(q);

%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%% compare different methods for calculating spkd (spike distance metric)
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------

% loop through units
for U = 1:1
	fprintf('********************************\n')
	fprintf('\tUNIT %d\n', U);
	fprintf('********************************\n')

	%-------------------------------------------------------------------------
	% extract the raw spiketrains for the current unit
	%-------------------------------------------------------------------------
	SpkTrains = AllSpkTrains((U*Conditions+1 - Conditions):(U*Conditions),:);
	
	%-------------------------------------------------------------------------
	% first, use original method from SpikeTrain2
	%-------------------------------------------------------------------------
	fprintf('\n%s\n', sepstr);
	fprintf('Original method:\n')
	fprintf('%s\n', sepstr);
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
		% loop through condition to use as "reference"
		for C1 = 1:Conditions
			% loop through trials
			for T1 = 1:Trials
				Dspike = zeros(Conditions, Trials);
				Dave = zeros(Conditions, 1);
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
				D(:, C1) = Dave';
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
	fprintf('%s\n', sepstr);
	fprintf('Original method, shuffled data:\n')
	fprintf('%s\n', sepstr);
	clear Ddistance Dspike Dave Dsum D Hvalue
	% preallocate things
	OrigShuf(U).Dconfuse = cell(Ncost, Nshuf);
	OrigShuf(U).Hvalue = zeros(Ncost, Nshuf);	
	% loop through cost values
	tic
	for Cindex = 1:Ncost
		Cost = q(Cindex);
		fprintf('q = %f\n', Cost);
		% loop through bootstrap count
		for B = 1:Nshuf
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
				end	% END T1 loop
			end	% END C1 loop
			% compute some stats from Ddistance (confusion) matrix
			OrigShuf(U).Dconfuse{Cindex, B} = Ddistance;
			OrigShuf(U).Hvalue(Cindex, B) = Hstat(Ddistance);
		end	% END Bootstrap loop
		fprintf('\tHvalue (mean): %f\n', mean(OrigShuf(U).Hvalue(Cindex, :)));
		fprintf('\tHvalue (std. dev): %f\n', std(OrigShuf(U).Hvalue(Cindex, :)));
	end	% END cost loop
	origbstraptime = toc;
	fprintf('original bootstrap took %f seconds\n', origbstraptime);
	fprintf('\n');

	
	%-------------------------------------------------------------------------
	% Use SPKDL MEX function instead of SPKD
	%-------------------------------------------------------------------------
	fprintf('\n%s\n', sepstr);
	fprintf('SPKDL method:\n')
	fprintf('%s\n', sepstr);
	clear Ddistance Dspike Dave Dsum D Hvalue
	% preallocate things
	MEXResults(U).Dconfuse = cell(Ncost, 1);
	MEXResults(U).Hvalue = zeros(Ncost, 1);	
	% loop through cost values	
	tic
	for Cindex = 1:Ncost
		Cost = q(Cindex);
		fprintf('q = %f\n', Cost);
		% allocate Distance matrix
		Ddistance = zeros(Conditions, Conditions);
		% loop through condition to use as "reference"
		for C1 = 1:Conditions
			% loop through trials
			for T1 = 1:Trials
				Dspike = zeros(Conditions, Trials);
				Dave = zeros(Conditions, 1);
				% for each pair of conditions (stimuli), compute spike
				% distance metric for all trials
				for C2 = 1:Conditions
					for T2=1:Trials
						Dspike(C2,T2) = (spkdmex(SpkTrains{C1,T1}', SpkTrains{C2,T2}', Cost)).^z;
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
				D(:, C1) = Dave';
			end;
		end;
		% compute some stats from Ddistance (confusion) matrix
		MEXResults(U).Dconfuse{Cindex} = Ddistance;
		MEXResults(U).Hvalue(Cindex) = Hstat(Ddistance);
		fprintf('\tHvalue: %f\n', MEXResults(U).Hvalue(Cindex));
		fprintf('\tDdistance = \n');
		disp(Ddistance);
	end	% END cost loop
	mex_time = toc;
	fprintf('spkdmex method took %f seconds\n', mex_time);
	fprintf('\n');
	
	%-------------------------------------------------------------------------
	% bootstrap the shuffled trains
	%-------------------------------------------------------------------------
	fprintf('%s\n', sepstr);
	fprintf('SPKD MEX method, shuffled data:\n')
	fprintf('%s\n', sepstr);
	clear Ddistance Dspike Dave Dsum D Hvalue
	% preallocate things
	MEXShuf(U).Dconfuse = cell(Ncost, Nshuf);
	MEXShuf(U).Hvalue = zeros(Ncost, Nshuf);	
	% loop through cost values	
	tic
	for Cindex = 1:Ncost
		Cost = q(Cindex);
		fprintf('q = %f\n', Cost);
		% loop through bootstrap count
		for B = 1:Nshuf
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
							Dspike(C2,T2) = (spkdmex(ShufTrains{C1,T1}, ShufTrains{C2,T2}, Cost)).^z;
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
				end	% END T1 loop
			end	% END C1 loop
			% compute some stats from Ddistance (confusion) matrix
			MEXShuf(U).Dconfuse{Cindex, B} = Ddistance;
			MEXShuf(U).Hvalue(Cindex, B) = Hstat(Ddistance);
		end	% END Bootstrap loop
		fprintf('\tHvalue (mean): %f\n', mean(MEXShuf(U).Hvalue(Cindex, :)));
		fprintf('\tHvalue (std. dev): %f\n', std(MEXShuf(U).Hvalue(Cindex, :)));
	end	% END cost loop
	mexb_time = toc;
	fprintf('bootstrap spkdmex method took %f seconds\n', mexb_time);
	fprintf('\n');
	fprintf('\n');	
	
	
	