close all; clear all; clear classes;

%------------------------------------------------------------
%------------------------------------------------------------
%% file names
%------------------------------------------------------------
%------------------------------------------------------------
matpath = '/Users/sshanbhag/Work/Data/DataWave/batmat'
matfilelist = {	'827_12-18-2012--2649_BBNrate_Sorted.mat', ...
						'829_01-10-2013--2859_strings block2_Sorted.mat', ...
						'827_12-18-2012--2954_syllable block_Sorted.mat', ...
					};
nMatfiles = length(matfilelist);
mIndx = 3;
matfile = fullfile(matpath, matfilelist{mIndx});
load(matfile)

%------------------------------------------------------------
%------------------------------------------------------------
%% options
%------------------------------------------------------------
%------------------------------------------------------------
probenum = 1;
unitnum = 255;
% window for psth.  time in milliseconds relative to stimulus onset timestamp
psthwin = [-100 1000];

%------------------------------------------------------------
%------------------------------------------------------------
%% create data object
%------------------------------------------------------------
%------------------------------------------------------------
d = DW.RateData(D, fullfile(matpath, matfile));

%------------------------------------------------------------
%------------------------------------------------------------
%% plot appropriate to type
%------------------------------------------------------------
%------------------------------------------------------------
% plot rasters/psth for each stimulus and attenuation
% individual stimuli will be plotted in separate figures (for all atten
% values)
d.plotRasterAndPSTH('probe', probenum, 'unit', unitnum, 'offset', [-100 0])
