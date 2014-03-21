clear all
close all
clear classes

%------------------------------------------------------------
%% file names
%------------------------------------------------------------
DLLName = 'C:\DataWave\DWShared\nsDWFile.dll';
inputpath = 'F:\Work\Data\DataWave\batmat\convertedDDF';
dobjoutpath = 'F:\Work\Data\DataWave\batmat\dobj';
% matfilelist = dir(fullfile(inputpath, '*.mat'));
matfilelist(1).name = '826_11-30-2012--3468_syllables block_Sorted.mat';
% matfilelist(1).name = '827_01-03-2013--2961_RepRate0_Sorted.mat';

inputpath = '/Users/sshanbhag/Work/Data/Bat/BatRestrainedData/NewMatFiles'
matfile = '07-25-2013-866-3172_StringsBlock.ddf_Sorted.mat'
matfile = '07-25-2013-866_2348_RepRate0_Sorted.mat'

%------------------------------------------------------------
%% create Data Object
%------------------------------------------------------------
n = 1;

% for n = 1:length(matfilelist)
	
	probenum = 1;
	unitnum = 255;
	psthwin = [-200 799];
	binsize = 10;
	%%
	
	matfile = fullfile(inputpath, matfile);
	[tmppath, tmpfile, tmpext] = fileparts(matfile);
	dobjfile = fullfile(dobjoutpath, [tmpfile '.mat']);

	% load matfile to get D 
	fprintf('%s: loading .mat file %s\n', mfilename, fullfile(inputpath, matfile));
	load(matfile);
	
	%% convert to obj
	d = DW.RateData(D, fullfile(inputpath, matfile));
	% plot unit waveforms (overlaid)
	d.plotUnitWaveforms('probe', probenum, 'unit', unitnum);

	%%
	clear PSTH bins spiketimes StimInfo   %in case old stored in Dobj file
	[PSTH, bins, spiketimes, StimInfo] = d.computePSTH(probenum, unitnum, binsize, psthwin);
   
	d.plotRasterAndPSTH('probe', probenum, 'unit', unitnum, 'offset', [-200 0])
	
	%------------------------------------------------------------
	%------------------------------------------------------------
	%% save Data Object
	%------------------------------------------------------------
	%------------------------------------------------------------
	% write the data object.  for debugging use, the D struct 
	% is also saved, but this will ultimately be redundant...
% 	save( fullfile(inputpath, ['Dobj_' matfile]), 'd', 'probenum', 'unitnum');


% end
