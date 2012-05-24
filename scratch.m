%-----------------------------------------------------------------------------
% working on object files
%-----------------------------------------------------------------------------

clear all

% create instance
d = DW.DWdata([pwd filesep 'BATtest.txt']);


[~, r, errflg] = d.readRawData;

d.parseMarkers(r);

% % set file
% d.file = [pwd filesep 'BATtest.txt'];
% 
% 
% [~, errFlg] = d.readDataWaveTextInfo;
% 
% [~, errFlg] = d.parseHeader
