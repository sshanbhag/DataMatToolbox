close all; clear all; clear classes;

%------------------------------------------------------------
%------------------------------------------------------------
%% file names
%------------------------------------------------------------
%------------------------------------------------------------
matpath = '/Users/sshanbhag/Work/Data/DataWave/batmat/convertedDDF';
matfilelist = {	'832_02-14-2013--2469_repRate20_Sorted.mat', ...
					};
nMatfiles = length(matfilelist);
mIndx = 1;
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
d = DW.ReprateData(D, fullfile(matpath, matfile));

%------------------------------------------------------------
%------------------------------------------------------------
%% plot appropriate to type
%------------------------------------------------------------
%------------------------------------------------------------
% plot rasters/psth for each stimulus and attenuation
% individual stimuli will be plotted in separate figures (for all atten
% values)
%d.plotRasterAndPSTH('probe', probenum, 'unit', unitnum, 'offset', [-100 0])

%% compute 1 ms bin psth
% [H, bins, spikes, stiminf] = d.computePSTH(1, 255, 10, [-100 900])

%% plot
% d.plotRasterAndPSTH('probe', probenum, 'unit', unitnum, 'offset', [-200 0])