clear all

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

%------------------------------------------------------------
%% init NS struct
%------------------------------------------------------------
[D, H] = DW.loadDDF(fullfile(datapath,filename));

%------------------------------------------------------------
%% get events
%------------------------------------------------------------
events = DW.readEvent(H, D);

%------------------------------------------------------------
%% get segments
%------------------------------------------------------------
segments = DW.readSegment(H, D);
return
%------------------------------------------------------------
%% get some analog data
%------------------------------------------------------------
analog = DW.readAnalog(H, D, [D.AnalogList(1) D.AnalogList(1)], [1 1001], [1000 1000]);

%------------------------------------------------------------
%% get neural entities (not fully implemented)
%------------------------------------------------------------
neural = DW.readNeural(H, D);
