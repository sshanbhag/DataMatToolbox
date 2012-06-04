%-----------------------------------------------------------------------------
% working on object files
%-----------------------------------------------------------------------------

clear all

% create instance
d = DW.DWdata([pwd filesep 'BATtest.txt']);

d.exportRawData

% load Markers - this will read in raw data from BATtest.txt and parse the
% marker information
d.loadData



