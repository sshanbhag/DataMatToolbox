%-----------------------------------------------------------------------------
% working on object files
%-----------------------------------------------------------------------------

clear all

% create instance
d = DW.DWdata([pwd filesep 'BATtest.txt']);

% read in raw text data.  Don't want to store this in the
% DWdata object itself
[~, r, errflg] = d.readRawData;

m = DW.Marker;


% d.parseMarkers(r);

% % set file
% d.file = [pwd filesep 'BATtest.txt'];
% 
% 
% [~, errFlg] = d.readDataWaveTextInfo;
% 
% [~, errFlg] = d.parseHeader
