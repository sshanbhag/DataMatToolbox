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


%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%% compare different methods for calculating spkd (spike distance metric)
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------

%-------------------------------------------------------------------------
% SDM settings
%-------------------------------------------------------------------------
% cost (q) in 1/seconds
q = 2.^-1;
% exponent
z = 1;

%-------------------------------------------------------------------------
% spike trains with distance = q (spiketimes in seconds)
%-------------------------------------------------------------------------
S1 = [0.5 1.0 5.0];
S2 = [0.5 (1.0 + q) 5.0];
S3 = [0.5 1.0 5.0 5.5];

%-------------------------------------------------------------------------
% settings for STAtoolkit
%-------------------------------------------------------------------------
opts.entropy_estimation_method = {'plugin'};
opts.shift_cost = q;
opts.start_time = 0;
opts.end_time = max(max([S1 S2])) + 10;
opts.metric_family = 0;
opts.clustering_exponent = z;	% normally this is -2. but for now, use 1
opts.parallel = 1;
opts.unoccupied_bins_strategy = -1;
Ncost = length(opts.shift_cost);

%-------------------------------------------------------------------------
% check distances using spkd
%-------------------------------------------------------------------------
% First, compute distance between identical spike trains (S1 and S1) using spkd
% (function from old toolkit (or from Kyle Nakamoto?)).  should be 0
d11.spkd = spkd(S1, S1, q);
% now with S1 and S2.  This should be equal to q, since the two trains only
% differ by an amount of q
d12.spkd = spkd(S1, S2, q);
% and with S1 and S3. should be 1, since there is one more spike in S3
d13.spkd = spkd(S1, S3, q);

%-------------------------------------------------------------------------
% now using STAtoolkit function, metric_distance
%-------------------------------------------------------------------------
% build Sin structure for STA toolkit - need to specify resolution, time, etc.
Sin1 = spikes2sta({S1; S1}, 'timescale', 1, 'timeresolution', 1, ...
							'start_time', opts.start_time, 'end_time', opts.end_time);
% compute metric, store results in tmp 
[tmp, tmpo] = metric(Sin1, opts);
% store tmp.d (distance matrix for the Spiketrains 1) in d__.sta
d11.sta = tmp.d;
% check if mismatch
if tmp.d(1, 2) ~= tmp.d(2, 1)
	warning('%s: mismatch in d matrix', mfilename);
	disp(d11.sta)
end

% same for S1 and S2
Sin2 = spikes2sta({S1; S2}, 'timescale', 1, 'timeresolution', 1, ...
							'start_time', opts.start_time, 'end_time', opts.end_time);
[tmp, tmpo] = metric(Sin2, opts);
d12.sta = tmp.d;
% check if mismatch
if tmp.d(1, 2) ~= tmp.d(2, 1)
	warning('%s: mismatch in d matrix', mfilename);
	disp(d12.sta)
end

% finally, S1, S3
Sin3 = spikes2sta({S1; S3}, 'timescale', 1, 'timeresolution', 1, ...
							'start_time', opts.start_time, 'end_time', opts.end_time);
[tmp, tmpo] = metric(Sin3, opts);
d13.sta = tmp.d;
% check if mismatch
if tmp.d(1, 2) ~= tmp.d(2, 1)
	warning('%s: mismatch in d matrix', mfilename);
	disp(d13.sta)
end

% display results
fprintf('\n\n');
fprintf('Spike train distances:\n')
fprintf('S1 -> S1:\n');
fprintf('\tspkd:\t%f\n', d11.spkd);
fprintf('\tsta:\t%f\n', d11.sta(1, 2));
fprintf('\n');

fprintf('S1 -> S2:\n');
fprintf('\tspkd:\t%f\n', d12.spkd);
fprintf('\tsta:\t%f\n', d12.sta(1, 2));
fprintf('\n');

fprintf('S1 -> S3:\n');
fprintf('\tspkd:\t%f\n', d13.spkd);
fprintf('\tsta:\t%f\n', d13.sta(1, 2));
fprintf('\n');


%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%% compare different methods for calculating Info for spkd 
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------

%-------------------------------------------------------------------------
% SDM settings
%-------------------------------------------------------------------------
% cost (q) in 1/seconds
q = 2.^-1;
% exponent
z = -2;

%-------------------------------------------------------------------------
% spike trains with distance = 0, q (spiketimes in seconds) and 2q
%-------------------------------------------------------------------------
S1 = [0.5 1.0 5.0];
S2 = [0.5 (1.0 + q) 5.0];
S3 = [0.5 (1.0 + q) (5.0 + q)];

%-------------------------------------------------------------------------
% first, use original method from SpikeTrain2
%-------------------------------------------------------------------------
Conditions = 3;
Trials = 50;
TotalUnits = 1;
% combine into SpkTrains
SpkTrains = cell(Conditions, Trials);
for c = 1:Conditions
	for t = 1:Trials
		switch c
			case 1
				SpkTrains{c, t} = S1;
			case 2
				SpkTrains{c, t} = S2;
			case 3
				SpkTrains{c, t} = S3;
		end
	end
end

fprintf('\n\n');
fprintf('Testing 3 trains/categories using original method\n')
% loop through cost values
for Cindex = 1:Ncost
	Cost = q(Cindex);
	fprintf('Cost = %f\n', Cost);
	% loop through units
	for Units = 1:TotalUnits
		% allocate Distance matrix
		Ddistance= zeros(Conditions, Conditions);
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
		Hvalue = Hstat(Ddistance);
		fprintf('Unit %d\n', Units)
		fprintf('\tHvalue: %f\n', Hvalue);
		fprintf('\tDdistance = \n');
		disp(Ddistance);		
	end;
end

%-------------------------------------------------------------------------
% then, use STA
%-------------------------------------------------------------------------
fprintf('Testing 3 trains/categories using STA toolkit\n')
% settings for STAtoolkit
opts.entropy_estimation_method = {'plugin'};
opts.shift_cost = q;
opts.start_time = 0;
opts.end_time = max(max([S1 S2])) + 10;
opts.metric_family = 0;
opts.clustering_exponent = z;
opts.parallel = 1;
opts.unoccupied_bins_strategy = -1;
Ncost = length(opts.shift_cost);

% build input and run
Sin = spikes2sta(SpkTrains, 'timescale', 1, 'timeresolution', 1, ...
							'start_time', opts.start_time, 'end_time', opts.end_time);
[tmp, tmpo] = metric(Sin, opts);

fprintf('H_sta: %f\n', tmp.table.information.value);
fprintf('Confusion Matrix:\n')
disp(tmp.cm);
fprintf('\n\n')



%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%% other synthetic data
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------

%-------------------------------------------------------------------------
% SDM settings
%-------------------------------------------------------------------------
% cost (q) in 1/seconds
q = 2.^-1;
% exponent
z = -2;

%-------------------------------------------------------------------------
% spike trains with distance = 0, q (spiketimes in seconds) and 2q
%-------------------------------------------------------------------------
S1 = [0.5 1.0 5.0];
S2 = [0.5 (1.0 + q) 5.0];
S3 = [0.5 (1.0 + q) (5.0 + q)];

%-------------------------------------------------------------------------
% first, use original method from SpikeTrain2
%-------------------------------------------------------------------------
Conditions = 3;
Trials = 50;
TotalUnits = 1;
% combine into SpkTrains
SpkTrains = cell(Conditions, Trials);
for c = 1:Conditions
	for t = 1:Trials
		r = 10*rand(1, 1);
		switch c
			case 1
				SpkTrains{c, t} = S1 + r;
			case 2
				SpkTrains{c, t} = S2 + r;
			case 3
				SpkTrains{c, t} = S3 + r;
		end
	end
end

fprintf('\n\n');
fprintf('Testing 3 trains (with random shift)/categories using original method\n')
clear Ddistance Dspike Dave Dsum D Hvalue
% loop through cost values
for Cindex = 1:Ncost
	Cost = q(Cindex);
	fprintf('Cost = %f\n', Cost);
	% loop through units
	for Units = 1:TotalUnits
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
		Hvalue = Hstat(Ddistance);
		fprintf('Unit %d\n', Units)
		fprintf('\tHvalue: %f\n', Hvalue);
		fprintf('\tDdistance = \n');
		disp(Ddistance);		
	end;
end

%-------------------------------------------------------------------------
% then, use STA
%-------------------------------------------------------------------------
fprintf('Testing 3 trains/categories using STA toolkit\n')
% settings for STAtoolkit
opts.entropy_estimation_method = {'plugin'};
opts.shift_cost = q;
opts.start_time = 0;
opts.end_time = max(max([S1 S2])) + 10;
opts.metric_family = 0;
opts.clustering_exponent = z;	
opts.parallel = 1;
opts.unoccupied_bins_strategy = -1;
Ncost = length(opts.shift_cost);

% build input and run
Sin = spikes2sta(SpkTrains, 'timescale', 1, 'timeresolution', 1, ...
							'start_time', opts.start_time, 'end_time', opts.end_time);
[tmp, tmpo] = metric(Sin, opts);

fprintf('H_sta: %f\n', tmp.table.information.value);
fprintf('Confusion Matrix:\n')
disp(tmp.cm);
fprintf('\n\n')

