close all; clear all;

%------------------------------------------------------------
%% file names
%------------------------------------------------------------
DLLName = 'C:\DataWave\DWShared\nsDWFile.dll';
datapath = 'F:\Work\Data\MG';
% filename = '12-12-2012--2854_BBNrate.ddf';
% filename = '12-12-2012--2854_FreqScan2.ddf';
% filename = '12-12-2012--2854_strings block.ddf';
% filename = '12-12-2012--2854_RepRate0.ddf';
filename = '01-03-2013--2961_syllable_block_new_sorted.ddf';
matfile = '01-03-2013--2961_syllable_block_new_sorted.mat';

%------------------------------------------------------------
%% convert to matfile or load matfile
%------------------------------------------------------------
if ~exist(fullfile(datapath, matfile), 'file')
	fprintf('%s: converting ddf file to mat file...\n', mfilename)
	D = DW.convertDDF2MAT(fullfile(datapath,filename), 'EVENT', 'SEGMENT', 'NEURAL');
	fprintf('%s: ....done\n', mfilename);
else
	fprintf('%s: loading .mat file %s\n', mfilename, fullfile(datapath, matfile));
	load(fullfile(datapath, matfile));
end

%------------------------------------------------------------
%% create Data
%------------------------------------------------------------
d = Data(D);