function Data = spikerater(infiles, spikewindows)
%------------------------------------------------------------------------
% spikerater.m
%------------------------------------------------------------------------
% script for spike rate computation and analysis from DataMat toolbox
%------------------------------------------------------------------------
%------------------------------------------------------------------------

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% define some constants/settings
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
Nwin = length(spikewindows);
Nfiles = length(infiles);


% pre-allocate Data cell array to hold data for each pair of BBN and LFH files
Data = cell(1, Nfiles);

% loop through file list
for fIndx = 1:Nfiles
	%------------------------------------------------------------------------
	% Load Data file
	%------------------------------------------------------------------------
	RAW = load(infiles{fIndx}, 'D', 'Stimulus', 'Background');

	fprintf('\n\n%d\n', fIndx);
	fprintf('\t%s\n',	RAW.D.Info.file);

	%--------------------------------------------------------------------------
	% now, need to parse filenames in order to extract condition
	% 		...BBN_..._1.mat  ---> pre odor (clean bedding)
	% 		...BBN_..._2.mat  ---> mild aversive odor (clean cotton)
	% 		...BBN_..._3.mat  ---> strong aversive odof (cat odor)
	%--------------------------------------------------------------------------
	% parse file name
	[~, tmpfname, ext] = fileparts(RAW.D.Info.file);
	% remove the _Sheetmaker from the filename
	tmpfname = regexprep(tmpfname, '(_Sheetmaker)', '');
	% remove the _spksorted from the filename
	tmpfname = regexprep(tmpfname, '(_spksorted)', '');
	% condition id value is last character in what remains of name
	RAW.D.Info.condition = str2num(tmpfname(end));

	fprintf('\tCondition: %d\n', RAW.D.Info.condition);

	% units match so set UnitList to one of the lists from the D structure
	Nunits = RAW.D.Info.Nunits;
	UnitList = RAW.D.Info.UnitList(:, 1:2);
	UnitIndex = RAW.D.Info.UnitList(:, 3);

	%--------------------------------------------------------------------------
	% determine number of Attenuation levels by looking at the Stimulus.Var
	% information in each condition's data structure
	%--------------------------------------------------------------------------
	% check number of attenuation values

	% # of attenuation values in RAW file
	nVar = length(RAW.Stimulus.Var);
	attenindex = 0;
	for n = 1:nVar
		% find R channel attenuation
		if strcmp(RAW.Stimulus.Var(n).name, 'AttenuationR')
			attenindex = n;
		end
	end
	% if found, store the values
	if attenindex
		attenVals = RAW.Stimulus.Var(attenindex).values;
	else
		% if not, error
		error('%s: attenuation not found in Var list for file %s', mfilename, RAW.D.Info.file);
	end
	Natten = length(attenVals);

	%------------------------------------------------------------------------
	% List of units to analyze
	%------------------------------------------------------------------------
	fprintf('\n# units found: %d\n', Nunits);
	fprintf('\t\tProbe\t\tCluster\n')
	for p = 1:Nunits
		fprintf('\t\t\t%d\t\t\t%d\n', UnitList(p, 1), UnitList(p, 2));
	end
	fprintf('Attenuation (dB):\t');
	fprintf('%d ', attenVals);	fprintf('\n');

	%--------------------------------------------------------------------------
	% construct data arrays for analysis
	%--------------------------------------------------------------------------
	% use buildSpikes function to build the spikes arrays
	%--------------------------------------------------------------------------
	RAW.Spikes = buildSpikes(RAW.Stimulus);

	%--------------------------------------------------------------------------
	%--------------------------------------------------------------------------
	% loop through units
	%--------------------------------------------------------------------------
	%--------------------------------------------------------------------------
	% clear UnitData variable
	clear UnitData;

	% initialize the unitIndex counter - this will be used to index the analyzed
	% data output in MUdata
	dIndex = 0;

	% loop through the units to compare for this file
	for u = 1:Nunits
		% increment counter
		dIndex = dIndex + 1;

		% get the unitNumber information
		unitNumber = UnitIndex(u, :);

		% store current unit information in UnitData structure
		UnitData(u).UnitNumber = unitNumber;
		UnitData(u).UnitInfo = RAW.D.Info.UnitInfo(unitNumber);
		UnitData(u).UnitList = UnitList(u, :);

		UnitData(u).spikes = cell(Natten, 1);
		UnitData(u).ntrials = zeros(Natten, 1);
		
		UnitData(u).Count_mean = cell(Natten, 1);
		UnitData(u).Count_std = cell(Natten, 1);
		UnitData(u).NetCounts = cell(Natten, 1);
		UnitData(u).Net_mean = zeros(Natten, 1);
		UnitData(u).Net_std = zeros(Natten, 1);
			
		for a = 1:Natten
			% get all spikes for current unit and attenuation
			UnitData(u).spikes{a} = RAW.Spikes{unitNumber, 1, a};
			% determine number of trials
			UnitData(u).ntrials(a) = length(RAW.Spikes{unitNumber, 1, a});

			% Assign total spike counts to the (surprise) NetSpikeCounts array

			% pre allocate SpikeCounts matrix
			NetSpikeCounts = zeros(1, UnitData(u).ntrials(a));
			for trial = 1:UnitData(u).ntrials(a)
				NetSpikeCounts(trial) = length(UnitData(u).spikes{a}{trial});
			end
			% store in overall data structure
			UnitData(u).NetCounts{a} = NetSpikeCounts;

			% compute mean and s.d. of net data
			UnitData(u).Net_mean(a) = mean(UnitData(u).NetCounts{a});
			UnitData(u).Net_std(a) = std(UnitData(u).NetCounts{a});

			%-------------------------------------------------------------------
			% get background spikes
			%-------------------------------------------------------------------
			UnitData(u).BG_timestamps = RAW.Background.Spiketimes{unitNumber(1)};
			UnitData(u).BG_count = length(RAW.Background.Spiketimes{unitNumber(1)});

			%-------------------------------------------------------------------
			%-------------------------------------------------------------------
			% Perform count of spikes in spike window
			%-------------------------------------------------------------------
			%-------------------------------------------------------------------
			tmpspikes = UnitData(u).spikes{a};
			%-------------------------------------------------------------------
			% Count number of spikes in spikewindows{}
			%-------------------------------------------------------------------
			for trial = 1:UnitData(u).ntrials(a)
				for win = 1:Nwin
					spikewin = spikewindows{win};
					% find spikes that are within the current window
					validSpikes = (spikewin(1) < tmpspikes{trial}) & ...
											(tmpspikes{trial} < spikewin(2));
					% store spiketimes
					if ~isempty(validSpikes)
						UnitData(u).spiketimes{a}{win, trial} = tmpspikes{trial}(validSpikes);
					else
						UnitData(u).spiketimes{a}{win, trial} = [];
					end
					% count them and store in count matrix
					UnitData(u).counts{a}(win, trial) = sum(validSpikes);
				end
			end

			UnitData(u).Count_mean{a} = zeros(Nwin, 1);
			UnitData(u).Count_std{a} = zeros(Nwin, 1);
			% compute mean and s.d.
			for w = 1:Nwin
				UnitData(u).Count_mean{a}(w) = mean(UnitData(u).counts{a}(w, :));
				UnitData(u).Count_std{a}(w) = std(UnitData(u).counts{a}(w, :));
			end

		end	% END OF a LOOP
	end	% END OF u LOOP
		
	%--------------------------------------------------------------------------
	%--------------------------------------------------------------------------
	% save data for this set of files
	%--------------------------------------------------------------------------
	%--------------------------------------------------------------------------
	clear tmp
	tmp.UnitIndex = UnitIndex;
	tmp.UnitData = UnitData;
	tmp.Info = RAW.D.Info;
	tmp.Stimulus = RAW.Stimulus;
	tmp.AttenVals = attenVals;

	Data{fIndx} = tmp;
	clear tmp;

end	% END OF fIndx LOOP




