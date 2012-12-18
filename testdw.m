
% data file name
% notes: 
%	04-23-2012_batstrings2.txt file is poorly-formed and will not be
%	read in correctly
%
%	2789string_full_LmarkerFilled_stringForWavefile.txt is from new
%	single-channel bat recordings
% 
%	BBNspikes.txt is broad band noise stimulus test file

fname = '2789string_full_LmarkerFilled_stringForWavefile.txt'

% clear old d object
clear d; 

% call constructor from DW package (DWdata) with filename, store object in d
d = DW.DWdata(fname)

% build stimulus structure
s = d.buildStimuliFromMarkers