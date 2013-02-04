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

%------------------------------------------------------------
%% init NS struct
%------------------------------------------------------------
ns = DW.NS(fullfile(datapath,filename));

%------------------------------------------------------------
%% get events
%------------------------------------------------------------
events = ns.getEvents;

%------------------------------------------------------------
%% get segments
%------------------------------------------------------------
for t = 1:1
	fprintf('t: %d\n', t);
	segments = ns.getSegment;
	timestamps = ns.getSegmentTimeStampsOnly;
end
return
%------------------------------------------------------------
%% get some analog data
%------------------------------------------------------------
analog = DW.readAnalog(H, D, [D.AnalogList(1) D.AnalogList(1)], [1 1001], [1000 1000]);

%------------------------------------------------------------
%% get neural entities (not fully implemented)
%------------------------------------------------------------
neural = DW.readNeural(H, D);
