function [out, errFlg] = parseDataWaveTextHeader(dwinfo)
%------------------------------------------------------------------------
% [out, errFlg] = parseDataWaveTextHeader(dwinfo)
%------------------------------------------------------------------------
% parse dwinfo header information
% 
% returns structure 
% 
%------------------------------------------------------------------------
% Input Arguments:
% 
% Output Arguments:
%
% 	out is a modified version of input dwinfo struct with following modifications
% 	and additions:
% 
% 	MarkerCols			indices of marker columns
% 	NMarkerCols			# of marker columns
% 	SpikeCols			indices of spike timestamp columns
% 	NSpikeCols			# of spike columns
% 	UnitCols				indices to unit id columns
% 	NUnitCols			# unit id columns (length(UnitCols))
% 	Ndatalines			# of data lines (adjusted for # of header lines)
% 							dwinfo.Nlines - N_HEADER_LINES
% 	MarkerTags			updated MarkerTags
% 
% 	errFlg	Error flag
% 					0		no error
% 					1		user cancelled file opening
% 					2		no fields found in header lines
% 					3		file not found
%
%------------------------------------------------------------------------
% See: readDataWaveHeader, loadDWfile
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 25 May, 2011 (SJS)
% 	- uses code snipped from readDataWaveTextInfo.m
%
% Revisions:
%	21 May, 2012 (SJS)
% 	 -	updated comments/documentation as move to obj. oriented design 
% 		progresses
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------


%-----------------------------------------------------------
% load defaults
%-----------------------------------------------------------
DataWaveDefaults;

%--------------------------------------------------------------------
% find spike timestamp header tags
%--------------------------------------------------------------------
% algorithm: search text fields in header using strncmp, then locate
%  non-zero elements in returned vector using find
%--------------------------------------------------------------------
timestamp_cols = strncmp(dwinfo.header.fields{1}, 'timestamp', length('timestamp'));
SpikeCols = find(timestamp_cols);
% make sure something was found
if isempty(SpikeCols)
	% if empty, return error
	errFlg = 2;
	out = [];
	return
end
NSpikeCols = length(SpikeCols);

%--------------------------------------------------------------------
% find unit header tags in header
%--------------------------------------------------------------------
tmp = strncmp(dwinfo.header.fields{1}, 'cluster', 1);
UnitCols = find(tmp);
% make sure something was found
if isempty(UnitCols)
	% if empty, warn user
	warning('DWFILE:TSTAMP', '%s: no unit/probe timestamp fields found in file %s header', ...
									mfilename, dwinfo.filename);
end
NUnitCols = length(UnitCols);


%--------------------------------------------------------------------
% find marker Timestamp header tags
%--------------------------------------------------------------------
% There are several ways to do this.  
% one that makes very few assumptions is to AND together the unit and spike
% timestamp cols, NOT them and those should be the marker cols.
%--------------------------------------------------------------------
tmp = strncmp(dwinfo.header.fields{1}, 'Timestamp', length('Timestamp'));
MarkerCols = find(tmp);
% make sure something was found
if isempty(MarkerCols)
	% if empty, warn user
	warning('DWFILE:MARKER', '%s: no Marker timestamp fields found in file %s header', ...
									mfilename, dwinfo.filename);
elseif length(MarkerCols) > 1
	% if unpredicted length, warn user
	warning('DWFILE:MARKER', '%s: %d Marker timestamp  fields found in file %s header', ...
									mfilename, length(MarkerCols), dwinfo.filename);
else
	% 3 of columns for marker information is the value of the first
	% spike timestamp column - 1 (all spike timestamps are written after
	% the marker information)
	NMarkerCols = SpikeCols(1) - 1;
end


%-----------------------------------------------------------
% get marker tags
%-----------------------------------------------------------
MarkerTags = dwinfo.header.fields{1}(1:NMarkerCols);

%-----------------------------------------------------------
% assign values to output structure
%-----------------------------------------------------------
out = dwinfo;
out.MarkerCols = MarkerCols;
out.NMarkerCols = NMarkerCols;
out.SpikeCols = SpikeCols;
out.NSpikeCols = NSpikeCols;
out.UnitCols = UnitCols;
out.NUnitCols = NUnitCols;
out.Ndatalines = dwinfo.Nlines - N_HEADER_LINES;
out.MarkerTags = MarkerTags;
errFlg = 0;
