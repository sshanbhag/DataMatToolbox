%------------------------------------------------------------------------
% spikerater.m
%------------------------------------------------------------------------
% script for spike rate computation and analysis from DataMat toolbox
%------------------------------------------------------------------------
% Some Notes:
%------------------------------------------------------------------------

clear all

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% define some options
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% size of time window after stimulus onset in which to count spikes
RESPONSEWINDOW_MS = 100;

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% define some constants/settings
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% default sound onset time time
soundOnsetTime = 0;


spikeCountWindow{1} = [0 RESPONSEWINDOW_MS];

Nwin = length(spikeCountWindow);


indata.filepath{1} = '/Users/sshanbhag/Work/Code/Matlab/dev/Analysis/DataWave/DataMatToolbox/2.2';
indata.filename{1} = '782_121511_17_3q_LFH_3_spksorted_Sheetmaker.mat';

%--------------------------------------------------------------------------
% define some plot constants/settings
%--------------------------------------------------------------------------

%------------------------------------------------------------------------
% List of files and units to analyze
%------------------------------------------------------------------------
%{
% open a dialog box to get data file name and path
[filename, pathname] = uigetfile('*sorted_Sheetmaker.mat', 'Select file for condition 1');
% if user clicks "Cancel", filename == 0, exit function
if ~filename
	return
end
% parse the returned path+filename into path, name and extension
[tmppath, tmpname, tmpext] = fileparts(fullfile(pathname, filename));
% build up the .MAT file names
indata.filesA = {[tmpname '.mat']};
% initialize path
indata.inputpathA = {tmppath};

% store current path
origpath = pwd;
% switch to data file A path
cd(tmppath);

% open a dialog box to get data file name and path
[filename, pathname] = uigetfile('*sorted_Sheetmaker.mat', 'Select file for condition 2');
% if user clicks "Cancel", filename == 0, exit function
if ~filename
	return
end
% parse the returned path+filename into path, name and extension
[tmppath, tmpname, tmpext] = fileparts(fullfile(pathname, filename));
% build up the .MAT file names
indata.filesB = {[tmpname '.mat']};
% initialize path
indata.inputpathB = {tmppath};

cd(origpath);

% load the mat files
try
	% load data from condition A
	A = load(fullfile(indata.inputpathA{1}, indata.filesA{1}), 'D', 'Stimulus');
catch
	error('%s: error in loading .mat file %s', indata.filesA{1}, fullfile(indata.inputpathA{1}, indata.filesA{1}));
end

try
	% load data from condition B
	B = load(fullfile(indata.inputpathB{1}, indata.filesB{1}), 'D', 'Stimulus');
catch
	error('%s: error in loading .mat file %s', indata.filesB{1}, fullfile(indata.inputpathB{1}, indata.filesB{1}));
end

%}


load(fullfile(indata.filepath{1}, indata.filename{1}), 'D', 'Stimulus', 'Background');

%------------------------------------------------------------------------
% List of units to analyze
%------------------------------------------------------------------------
disp(sprintf('File %s:\t\t %d units found', indata.filename{1}, D.Info.Nunits))


% type code for units:
%	
% 	Value		String
% 	 1			 junk
% 	 2			 messy
% 	 3			 multi
% 	 4			 sorted
%	 5			 unclassified
type_strings = { ...
	'junk'			, ...
	'messy'			, ...
	'multi'			, ...
	'sorted'			, ...
	'unclassified'	 ...
};

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
	
Nunits = D.Info.Nunits;
	
%--------------------------------------------------------------------------
% determine number of Attenuation levels by looking at the Stimulus.Var
% information in each condition's data structure
%--------------------------------------------------------------------------

% check number of attenuation values
nVar = length(Stimulus.Var);
attenindex = 0;
for n = 1:nVar
	if strcmp(Stimulus.Var(n).name, 'AttenuationR')
		attenindex = n;
	end
end
if attenindex
	indata.attenVals = Stimulus.Var(attenindex).values;
else
	error('%s: attenuation not found in Var list for file %s', mfilename, indata.filename{1});
end
Natten = length(indata.attenVals);
	
%--------------------------------------------------------------------------
% construct data arrays for mutual information analysis
%--------------------------------------------------------------------------
% use buildSpikes function to build the spikes arrays
%--------------------------------------------------------------------------
Spikes = buildSpikes(Stimulus);
	
%--------------------------------------------------------------------------
% determine number of trials (use second atten value to eliminate extra
% trial at end)
%--------------------------------------------------------------------------
for n = 1:Natten
	indata.ntrials(n) = length(Spikes{1, 1, n});
end

NtrialList = indata.ntrials;
	

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% loop through units
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% clear UnitData variable
clear UnitData;
	
% initialize the unitIndex counter - this will be used to index the analyzed
% data output in MUdata
unitIndex = 0;

% loop through the units to compare for this pair of recordings (files)
for unitNumber = 1:Nunits
	% increment counter
	unitIndex = unitIndex + 1;
		
	% store current unit information in UnitData structure
	UnitData(unitIndex).UnitNumber = unitNumber;
	UnitData(unitIndex).UnitInfo = D.Info.UnitInfo(unitNumber);
	
	%--------------------------------------------------------------------------
	%--------------------------------------------------------------------------
	% Loop through attenuation values
	%--------------------------------------------------------------------------
	%--------------------------------------------------------------------------
	for attIndex = 1:Natten
		disp(sprintf('...Atten = %d dB', indata.attenVals(attIndex)))
			
		% get all spikes for current unit and attenuation, condition A
		unitresp = Spikes{unitNumber, 1, attIndex};
		UnitData(unitIndex).unitresp{attIndex} = unitresp;
		
		% Assign spike counts to the (surprise) SpikeCounts array
		% assuming that both conditions have same number of trials... 
	 	% pre allocate SpikeCounts matrix
	 	NetSpikeCounts = zeros(1, NtrialList(attIndex));
		for trial = 1:NtrialList(attIndex)
			NetSpikeCounts(1, trial) = length(unitresp{trial});
		end
			
		% store in overall data structure
		UnitData(unitIndex).NetSpikeCounts{attIndex} = NetSpikeCounts;

		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		% Count number of spikes in spikeCountWindow{}
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		for trial = 1:NtrialList(attIndex)
			for win = 1:Nwin
				spikewin = spikeCountWindow{win};
				% find spikes that are within the current window
				validSpikes = (spikewin(1) < unitresp{trial}) & ...
									(unitresp{trial} < spikewin(2));
									
				% count them and store in count matrix
				UnitData(unitIndex).counts{win}(attIndex, trial) = sum(validSpikes);
			end
		end
						
	end		% END OF attenIndex LOOP
	
	
	%------------------------------------------------------------------------
	% get background spikes
	%------------------------------------------------------------------------
	
	UnitData(unitIndex).BG_timestamps = Background.Spiketimes{unitIndex};
	UnitData(unitIndex).BG_count = length(Background.Spiketimes{unitIndex});
	
end		% END OF unitIndex LOOP



%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% process data
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%{
for indexF = 1:Nfiles
	UnitData = UnitDataList{indexF};

	for u = 1:length(indata.units_to_compare{indexF})

		for a = 1:Natten

			% loop through the different spike rate windows
			for w = 1:length(spikeCountWindow)
				% get mean and stddeviation of spike counts
				UnitData(u).meanCountsA{w, a} = mean(UnitData(u).countsA{w}(a, :));
				UnitData(u).stddevCountsA{w, a} = std(UnitData(u).countsA{w}(a, :));

				UnitData(u).meanCountsB{w, a} = mean(UnitData(u).countsB{w}(a, :));
				UnitData(u).stddevCountsB{w, a} = std(UnitData(u).countsB{w}(a, :));
			end
			
			% find peak information
			[maxval maxindex] = max(UnitData(u).I_sliding{a});
			UnitData(u).peakMI(a) = maxval;
			UnitData(u).peakMI_time(a) = maxindex * windowDuration;
			
		end
		
		% find peak MI across all attenuation values
		[maxval maxindex] = max(UnitData(u).peakMI);
		UnitData(u).globalPeakMI = maxval;
		UnitData(u).globalPeakMI_Time = UnitData(u).peakMI_time(maxindex);
		UnitData(u).globalPeakMI_Atten = indata.attenVals{indexF}(maxindex);
	end
	
	MIData{indexF} = UnitData;
end
%}

for u = 1:Nunits
	for a = 1:Natten
		UnitData(u).Net_mean(a) = mean(UnitData(u).NetSpikeCounts{a});
		UnitData(u).Net_std(a) = std(UnitData(u).NetSpikeCounts{a});
		
		
		for w = 1:length(spikeCountWindow)
			UnitData(u).Count_mean{w, a} = mean(UnitData(u).counts{w}(a, :));
			UnitData(u).Count_std{w, a} = std(UnitData(u).counts{w}(a, :));
		end
		
	end
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% export data
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%----------------------------------------------------------------
% Write all data to MutInf mat file and MI data to _MI.txt file
%----------------------------------------------------------------
%{
for indexF = 1:Nfiles
	attenVals = indata.attenVals{indexF};
	Natten = indata.Natten{indexF};
	
	endA = strfind(indata.filesA{indexF}, 'sorted') - 2;
	endB = strfind(indata.filesB{indexF}, 'sorted') - 2;
	
	outname = [indata.filesA{indexF}(1:endA) '-' indata.filesB{indexF}(1:endB) '_MI.txt'];
	
	matfilename = [indata.filesA{indexF}(1:endA) '-' indata.filesB{indexF}(1:endB) '_MutInf.mat'];
	disp(['Writing all data to file ' matfilename '...']);
	save(matfilename, 'MIData', '-MAT')
	disp('...done')
	
	disp(['Exporting MI data to text file ' outname '...']);
	fp = fopen(outname, 'w');
	
	% write headers
	fprintf(fp, 'FileA\t');
	fprintf(fp, 'FileB\t');
	fprintf(fp, 'unitNumber\t');
	fprintf(fp, 'unitType\t');
	fprintf(fp, 'unitTypeString\t');
	% I_global for all attenuation values
	for a = 1:Natten
		astr = sprintf('I_global(%d dB)', attenVals(a));
		fprintf(fp, '%s\t', astr);
	end
	% I_sh for all attenuation values
	for a = 1:Natten
		astr = sprintf('I_sh(%d dB)', attenVals(a));
		fprintf(fp, '%s\t', astr);
	end
	% peakMI for all attenuation values
	for a = 1:Natten
		astr = sprintf('peakMI(%d dB)', attenVals(a));
		fprintf(fp, '%s\t', astr);
	end
	% peakMI_sh for all attenuation values
	for a = 1:Natten
		astr = sprintf('peakMI_time(%d dB)', attenVals(a));
		fprintf(fp, '%s\t', astr);
	end
	% terminate header line
	fprintf(fp, '\n');
	
	% write data for all units
	UnitData = MIData{indexF};
	% loop through units
	for u = 1:length(UnitData)
		Utmp = UnitData(u);
		fprintf(fp, '%s\t', indata.filesA{indexF});
		fprintf(fp, '%s\t', indata.filesB{indexF});
		fprintf(fp, '%d\t', Utmp.unitNumber);
		fprintf(fp, '%d\t', Utmp.unitType);
		fprintf(fp, '%s\t', Utmp.unitTypeString);
		% I_global for all attenuations
		for a = 1:Natten
			fprintf(fp, '%f\t', Utmp.I_global(a));
		end
		% I_sh for all attenuation values
		for a = 1:Natten
			fprintf(fp, '%f\t', Utmp.Ish_global(a));
		end
		% peakMI for all attenuations
		for a = 1:Natten
			fprintf(fp, '%f\t', Utmp.peakMI(a));
		end
		% peakMI time for all attenuations
		for a = 1:Natten
			fprintf(fp, '%f\t', Utmp.peakMI_time(a));
		end
		% end of line
		fprintf(fp, '\n');
	end		% END OF u LOOP


	% close file
	fclose(fp);
	disp('...done')

end		% END of indexF LOOP

%}

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% write to rate file
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%{
for indexF = 1:Nfiles
	UnitData = MIData{indexF};

	attenVals = indata.attenVals{indexF};
	Natten = indata.Natten{indexF};
	
	endA = strfind(indata.filesA{indexF}, 'sorted') - 2;
	endB = strfind(indata.filesB{indexF}, 'sorted') - 2;
	
	outname = [indata.filesA{indexF}(1:endA) '-' indata.filesB{indexF}(1:endB) '_RATE.txt'];
	disp(['Exporting SpikeCount data to text file ' outname '...'])
	fp = fopen(outname, 'w');
	
	% write headers
	fprintf(fp, 'FileA\t');
	fprintf(fp, 'FileB\t');
	fprintf(fp, 'unitNumber\t');
	fprintf(fp, 'unitType\t');
	fprintf(fp, 'unitTypeString\t');
	
	% loop through spike count windows
	for w = 1:length(spikeCountWindow)
		% meanCountsA for all attenuation values
		for a = 1:Natten
			astr = sprintf('meanCountsA(%d db %d-%d ms)', ...
									attenVals(a), spikeCountWindow{w}(1), spikeCountWindow{w}(2));
			fprintf(fp, '%s\t', astr);
		end
		% stddevCountsA for all attenuation values
		for a = 1:Natten
			astr = sprintf('stddevCountsA(%d db %d-%d ms)', ...
									attenVals(a), spikeCountWindow{w}(1), spikeCountWindow{w}(2));
			fprintf(fp, '%s\t', astr);
		end
		% meanCountsB for all attenuation values
		for a = 1:Natten
			astr = sprintf('meanCountsB(%d db %d-%d ms)', ...
									attenVals(a), spikeCountWindow{w}(1), spikeCountWindow{w}(2));
			fprintf(fp, '%s\t', astr);
		end
		% stddevCountsB for all attenuation values
		for a = 1:Natten
			astr = sprintf('stddevCountsB(%d db %d-%d ms)', ...
									attenVals(a), spikeCountWindow{w}(1), spikeCountWindow{w}(2));
			fprintf(fp, '%s\t', astr);
		end
	end
	% end of header line
	fprintf(fp, '\n');
	
	% write data
	for u = 1:length(UnitData)
		Utmp = UnitData(u);
		fprintf(fp, '%s\t', indata.filesA{indexF});
		fprintf(fp, '%s\t', indata.filesB{indexF});
		fprintf(fp, '%d\t', Utmp.unitNumber);
		fprintf(fp, '%d\t', Utmp.unitType);
		fprintf(fp, '%s\t', Utmp.unitTypeString);
		
		% loop through spike count windows
		for w = 1:length(spikeCountWindow)

			for a = 1:Natten
				fprintf(fp, '%f\t', Utmp.meanCountsA{w, a});
			end

			for a = 1:Natten
				fprintf(fp, '%f\t', Utmp.stddevCountsA{w, a});
			end

			for a = 1:Natten
				fprintf(fp, '%f\t', Utmp.meanCountsB{w, a});
			end

			for a = 1:Natten
				fprintf(fp, '%f\t', Utmp.stddevCountsB{w, a});
			end
		end
		% end of line
		fprintf(fp, '\n');
	end
	fclose(fp);
	disp('...done')
end		%END OF FILE LOOP

%}


fp = fopen('tmpout.csv', 'w');

fprintf(fp, 'SpikeCountWindowLength,%f,msec,\n', RESPONSEWINDOW_MS);

fprintf(fp, 'File,');
fprintf(fp, 'UnitNum,');
fprintf(fp, 'Probe,');
fprintf(fp, 'Cluster,');
fprintf(fp, 'BGcount,');
for a = 1:Natten
	fprintf(fp, '%d dB Net_mean,', indata.attenVals(a));
end
for a = 1:Natten
	fprintf(fp, '%d dB Net_sd,', indata.attenVals(a));
end
for a = 1:Natten
	fprintf(fp, '%d dB Count_mean,', indata.attenVals(a));
end
for a = 1:Natten
	fprintf(fp, '%d dB Count_sd,', indata.attenVals(a));
end

fprintf(fp, '\n');

for u = 1:Nunits
	fprintf(fp, '%s,', indata.filename{1});
	fprintf(fp, '%d,', UnitData(u).UnitInfo.unit);
	fprintf(fp, '%d,', UnitData(u).UnitInfo.probe);
	fprintf(fp, '%d,', UnitData(u).UnitInfo.cluster);
	fprintf(fp, '%f,', UnitData(u).BG_count);
	
	for a = 1:Natten
		fprintf(fp, '%f,', UnitData(u).Net_mean(a));
	end
	for a = 1:Natten
		fprintf(fp, '%f,', UnitData(u).Net_std(a));
	end
	for a = 1:Natten
		fprintf(fp, '%f,', UnitData(u).Count_mean{1, a});
	end
	for a = 1:Natten
		fprintf(fp, '%f,', UnitData(u).Count_std{1, a});
	end
	
	
	
	fprintf(fp, '\n');
end

if fp ~= 1
	fclose(fp);
end

%----------------------------------------------------------------
%----------------------------------------------------------------
% Write peak MI data
%----------------------------------------------------------------
%----------------------------------------------------------------
%{
for indexF = 1:Nfiles
	endA = strfind(indata.filesA{indexF}, 'sorted') - 2;
	endB = strfind(indata.filesB{indexF}, 'sorted') - 2;
	
	outname = [indata.filesA{indexF}(1:endA) '-' indata.filesB{indexF}(1:endB) '_GlobalPeakMI.txt'];
	

	disp(['Exporting Global Peak MI data to text file ' outname '...']);
	fp = fopen(outname, 'w');
	
	% write headers
	fprintf(fp, 'FileA\t');
	fprintf(fp, 'FileB\t');
	fprintf(fp, 'unitNumber\t');
	fprintf(fp, 'unitType\t');
	fprintf(fp, 'unitTypeString\t');
	% globalPeakMI 
	fprintf(fp, 'globalPeakMI(bits)\t');
	% globalPeakMI_Time
	fprintf(fp, 'globalPeakMI_Time(ms)\t');
	% globalPeakMI_Attenuation
	fprintf(fp, 'globalPeakMI_Atten(dB)\t');
	% end of header line
	fprintf(fp, '\n');
	
	% write data for all units
	UnitData = MIData{indexF};
	% loop through units
	for u = 1:length(UnitData)
		Utmp = UnitData(u);
		fprintf(fp, '%s\t', indata.filesA{indexF});
		fprintf(fp, '%s\t', indata.filesB{indexF});
		fprintf(fp, '%d\t', Utmp.unitNumber);
		fprintf(fp, '%d\t', Utmp.unitType);
		fprintf(fp, '%s\t', Utmp.unitTypeString);
		fprintf(fp, '%f\t', UnitData(u).globalPeakMI);
		fprintf(fp, '%f\t', UnitData(u).globalPeakMI_Time);
		fprintf(fp, '%f\t', UnitData(u).globalPeakMI_Atten);
		% end of line
		fprintf(fp, '\n');
	end		% END OF u LOOP
	% close file
	fclose(fp);
	disp('...done')

end		% END of indexF LOOP
%}
