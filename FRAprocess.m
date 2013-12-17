%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------
% FRAprocess.m
%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------

%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------
%% initial setup
%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------

% clear all vars and close all plots to avoid confusion
% close all; 
% clear all;
FORCEFLAG = 0;
SAVEFLAG = 0;
SPIKE_WINDOW = [0 250];
FULL_WINDOW = [0 1000];

%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------
%% file names
%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------
%
%-----------------------------------------------------------------------------
%	1	'829_01-05-2013--2729_FRA_Sorted.mat'
% 	2	'829_01-07-2013--2535_FRA_Sorted.mat'
% 	3	'829_01-07-2013--2388_FRA_Sorted.mat'
% 	4	'829_01-09-2013--2893_FRA_Sorted.mat'
% 	5	'829_01-10-2013--2680_FRA_Sorted.mat'
% 	6	'829_01-14-2013--2579_FRA_Sorted.mat'
% 	7	'829_01-14-2013--2579_FRA_Sorted.mat'
% 	8	'832_01-25-2013--2714_FRA_Sorted.mat'
% 	9	'832_01-25-2013--2774_FRA_Sorted.mat'
% 	10	'832_01-28-2013--2864_FRA_Sorted.mat'
% 	11	'835_02-05-2013--3061_FRA_Sorted.mat'
% 	12	'835_02-05-2013--3061_FRA_Sorted.mat'
% 	13	'832_02-13-2013--3055_FRA_Sorted.mat'
% 	14	'835_02-28-2013--3042_FRA_Sorted.mat'
%-----------------------------------------------------------------------------

%------------------------------------------------------------
% set path and library depending on system
%------------------------------------------------------------
switch computer
	case {'PCWIN', 'PCWIN64'}
		DLLName = 'C:\DataWave\DWShared\nsDWFile.dll';
		datapath = 'F:\Work\Data\DataWave\DDFtest';
	case 'MACI64'
		datapath = '/Users/sshanbhag/Work/Data/Bat/BatRestrainedData/FRAdata';
end

addpath('/Users/sshanbhag/Work/Code/Matlab/dev/Toolboxes/SpikeUtilities/2.0');

%------------------------------------------------------------
% information file
%------------------------------------------------------------
infofile = 'FRA_UnitInformation.csv';

% load information
fp = fopen(fullfile(datapath, infofile), 'r');
tmp = textscan(fgetl(fp), '%s', 'Delimiter', ',');
fieldnames = tmp{1};
nfields = length(fieldnames);
FileData = cell(1, nfields);
lIndex = 1;
while(~feof(fp))
	tmp = textscan(fgetl(fp), '%s', 'Delimiter', ',');
	tmpl = tmp{1}';
	% append empty cells if necessary
	if length(tmpl) < nfields
		tmpl = [tmpl cell(1, nfields - length(tmpl))];
	end
	FileData(lIndex, :) = tmpl;
	lIndex = lIndex + 1;
end
fclose(fp);
FRAinfo = cell2struct(FileData', fieldnames);
clear tmp tmpl lIndex
nFiles = length(FRAinfo);

%------------------------------------------------------------
%% Process and plot data files
%------------------------------------------------------------
% freq/response area, sorted spikes (bat)
% for n = 1:nFiles
for n = 8:8
	clear D d S basename matfile figfile pngfile fra
	
	fra = FRAinfo(n);
	basename = [fra.Animal '_' fra.Date '--' fra.Depth]
	matfile = [basename '_FRA_Sorted.mat']
	figfile = [basename '_FRA.fig'];
	pngfile = [basename '_FRA.png'];

	%-----------------------------------------------------------------------------
	%-----------------------------------------------------------------------------
	%% convert DDF data to matfile (via NeuroShare) OR load matfile
	%-----------------------------------------------------------------------------
	%-----------------------------------------------------------------------------
	if ~exist(fullfile(datapath, matfile), 'file') || FORCEFLAG
		fprintf('%s: converting ddf file to mat file...\n', mfilename)
		try
			D = DW.convertDDF2MAT(fullfile(datapath,filename), 'EVENT', 'SEGMENT', 'NEURAL');
		catch errMsg
			errMsg %#ok<NOPTS>
		end
		fprintf('%s: ....done\n', mfilename);
	else
		fprintf('%s: loading .mat file %s\n', mfilename, fullfile(datapath, matfile));
		load(fullfile(datapath, matfile));
	end

	%-----------------------------------------------------------------------------
	%-----------------------------------------------------------------------------
	%% create Data Object from DDF neuroshare data matfile
	%-----------------------------------------------------------------------------
	%-----------------------------------------------------------------------------
	d = DW.FRAdata(D, fullfile(datapath, matfile));

	%-----------------------------------------------------------------------------
	%-----------------------------------------------------------------------------
	%% plot FRA data
	%-----------------------------------------------------------------------------
	%-----------------------------------------------------------------------------
	d.plotFRA(1, 255, SPIKE_WINDOW)
	fprintf('%s has %d sweeps per stimulus\n\n', matfile, d.Stimuli.sweepsList);
	
	%-----------------------------------------------------------------------------
	%-----------------------------------------------------------------------------
	%% plot PSTH
	%-----------------------------------------------------------------------------
	%-----------------------------------------------------------------------------
% 	S = d.getFRAspikes('probe', 1, 'unit', 255, 'window', SPIKE_WINDOW);
	% get the spiketimes sorted by attenuation and frequency
	sortS = d.getSortedFRASpikes(1, 255, FULL_WINDOW);
	% store number of sweeps per atten/freq combination
	nsweeps(n) = d.Stimuli.sweepsList;
	% allocate storage for data
	%	spikes_ms	{Natten, Nfreqs}{nsweeps} array of spike times, milliseconds
	%	spikes_s		{Natten, Nfreqs}{nsweeps} spike times, seconds
	%	spikes_sdf	{Natten, Nfreqs}{nsweeps} spike density function
	%	spikes_sdfsum	{Natten, Nfreqs} summed spike density function
	spikes_ms = cell(d.Natten, d.Nfreqs);
	spikes_s = cell(d.Natten, d.Nfreqs);
	spikes_sdf = cell(d.Natten, d.Nfreqs);
	spikes_sdfsum = cell(d.Natten, d.Nfreqs);
	FRApsth = cell(d.Natten, d.Nfreqs);
	FRApsth_bins = [];
	
	psthbinwidth = 20;
	psthtime = 1000;
	
	sdf_duration = 1;	% max duration of sdf in seconds
	sdf_Fs = 2000;
	prise = ms2samples(2, sdf_Fs);
	pdecay = ms2samples(5, sdf_Fs);

	% convert spike times in microseconds to spiketimes in msec and seconds
	for f = 1:d.Nfreqs
		for a = 1:d.Natten
			spikes_ms{a, f} = cell(nsweeps(n), 1);
			spikes_s{a, f} = cell(nsweeps(n), 1);
			spikes_sdf{a, f} = cell(nsweeps(n), 1);
			for s = 1:nsweeps(n)
				spikes_ms{a, f}{s} = 0.001 * sortS.SpikeTimes{a, f}{s};
				spikes_s{a, f}{s} = 0.001 * spikes_ms{a, f}{s};
				spikes_sdf{a, f}{s} = poissconv(int32(sdf_Fs * spikes_s{a, f}{s}), ...
															prise, pdecay, sdf_Fs * sdf_duration);
			end
			% compute average sdf
 				spikes_sdfsum{a, f} = sum(cell2mat(spikes_sdf{a, f}));
			if isempty(FRApsth_bins)
				[FRApsth{a, f}, FRApsth_bins] = ...
											psth(spikes_ms{a, f}, psthbinwidth, psthtime);
			else
				FRApsth{a, f} = psth(spikes_ms{a, f}, psthbinwidth, psthtime);
			end
		end	% END a LOOP
	end	% END f LOOP
	
	figure(2)
	subplot(211)
	[h, bins] = psth(spikes_ms{1, 1}, psthbinwidth, psthtime);
	bar(bins, h)
	subplot(212)
	bar(FRApsth_bins, FRApsth{1, 1})
	
	%% plot sdfs
	figure(3)
	tvec = ( (1:length(spikes_sdfsum{1, 1})) - 1 ) ./ sdf_Fs;
	pindex = 1;
	for a = 1:d.Natten
		for f = 1:d.Nfreqs
			subplot(d.Natten, d.Nfreqs, pindex)
			plot(tvec, spikes_sdfsum{a, f});
			xlim([0 1]);
			ylim([0 5]);
			box off
			
			if (f ~= 1) || (a ~= d.Natten)
				set(gca, 'XTickLabel', '');
				set(gca, 'YTickLabel', '');
			end
			set(gca, 'TickDir', 'out')
			pindex = pindex + 1;
		end	% END f LOOP
	end	% END a LOOP
	set(gcf, 'Color', get(gca, 'Color'));
	
	%-----------------------------------------------------------------------------
	%-----------------------------------------------------------------------------
	%% plot FRA data
	%-----------------------------------------------------------------------------
	%-----------------------------------------------------------------------------
	d.plotFRA(1, 255, SPIKE_WINDOW)
	fprintf('%s has %d sweeps per stimulus\n\n', matfile, d.Stimuli.sweepsList);
	
	%-----------------------------------------------------------------------------
	%-----------------------------------------------------------------------------
	%% save Data Object
	%-----------------------------------------------------------------------------
	%-----------------------------------------------------------------------------
	if SAVEFLAG
		save(fullfile(datapath, ['Dobj_' matfile]), 'd');
		saveas(gcf, fullfile(datapath, figfile), 'fig');
		saveas(gcf, fullfile(datapath, pngfile), 'png');
	end
end