clear all
close all

DataWaveDefaults;

%------------------------------------------------------------
%% file names
%------------------------------------------------------------
DLLName = 'C:\DataWave\DWShared\nsDWFile.dll';
datapath = 'F:\Work\Data\MG';
filename = '12-12-2012--2854_BBNrate.ddf';
% filename = '12-12-2012--2854_FreqScan2.ddf';
% filename = '12-12-2012--2854_strings block.ddf';
% filename = '12-12-2012--2854_RepRate0.ddf';
% filename = '01-03-2013--2961_syllable_block_new_sorted.ddf';

%------------------------------------------------------------
%% init DW struct
%------------------------------------------------------------
D = DW.DWdata(fullfile(datapath,filename));

%------------------------------------------------------------
%% load markers
%------------------------------------------------------------
[~, Events] = D.loadMarkers;

%------------------------------------------------------------
%% load stimuli
%------------------------------------------------------------
D.loadStimuli;