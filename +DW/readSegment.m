function [d, s] = readSegment(Hfile, Data, varargin)
%------------------------------------------------------------------------
%[d, s] = readSegment(Hfile, Data)
%[d, s] = readSegment(Hfile, Data, 'Segments', SegmentList)
%[d, s] = readSegment(Hfile, Data, 'LoadWaves', <<1== yes, 2 == no>>)
%------------------------------------------------------------------------
% 	Segment Entities (3)
% 		Short, time-stamped segments of digitized analog
% 		signals in which the segments are separated by variable amounts of time.
% 		Segment Entities can contain data from more than one source. They are
% 		intended to represent discontinuous analog signals such as extracellular
% 		spike waveforms from electrodes or groups of electrodes.
% 		Each index of a segment entity refers to a short, time-stamped
% 		segment of analog data from one or more sources.  The number of
% 		indexes is equal to the number of entries for that segment entity in
% 		the file.
% 
%------------------------------------------------------------------------
% Input Arguments:
%
% Output Arguments:
%------------------------------------------------------------------------
% See also: Neuroshare MATLAB API
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 2 January, 2013 (SJS)
%
% Revisions:
%	10 Jan 2013 (SJS): renamed to readSegment
%	29 Jan 2013 (SJS): 
% 	 -	added to +DW package
% 	 -	added option for SegmentList to be provided as input
% 
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

% check if there are any segment entities in the Data structure
if Data.nSegment == 0
	disp('No segment entities available');
	d = [];
	return
end

% defaults
SegmentList = Data.SegmentList;
LoadWaves = 1;

% set from options
argn = 1;
while argn <= length(varargin)
	switch upper(varargin{argn})
		case 'SEGMENTS'
			SegmentList = varargin{argn+1};
			argn = argn + 2;
		case 'LOADWAVES'
			LoadWaves = varargin{argn+1};
			argn = argn + 2;
		otherwise
			error('%s: unknown option %s', mfilename, varargin{argn});
	end		
end

% check # of segments to load
nSegments = length(SegmentList);
if ~nSegments
	fprintf('ERROR(%s): no segments in SegmentList\n\n', mfilename);
	d = [];
	return
end

% loop through segments
s = zeros(nSegments, 1);
for n = 1:nSegments
	% get info
	[s(n), tmp.Info] = ns_GetSegmentInfo(Hfile, SegmentList(n));
	[s(n), tmp.SourceInfo] = ns_GetSegmentSourceInfo(Hfile, SegmentList(n), 1);
	% get # of items (waveforms) for this segment
	tmp.ItemCount = Data.EntityInfo(SegmentList(n)).ItemCount;
	if isempty(tmp.ItemCount) || (tmp.ItemCount == 0)
		tmp.TimeStamp = [];
		tmp.WaveForm = {};
		tmp.Nsamples = [];
		tmp.UnitID = [];
	else
		% pre-allocate storage
		tmp.TimeStamp = zeros(tmp.ItemCount, 1);
		if LoadWaves
			tmp.WaveForm = cell(tmp.ItemCount, 1);
			tmp.Nsamples = zeros(tmp.ItemCount, 1);
		else
			tmp.WaveForm = {};
			tmp.Nsamples = {};
		end
		tmp.UnitID = zeros(tmp.ItemCount, 1);

		% two different options depending on LoadWaves value
		if LoadWaves
			for m = 1:tmp.ItemCount
				% Load all data (including waveforms) for selected channel
				[s(n), tmp.TimeStamp(m), tmp.WaveForm{m}, tmp.Nsamples(m), tmp.UnitID(m)] = ...
						ns_GetSegmentData(Hfile, SegmentList(n), m);
			end
		else
			for m = 1:tmp.ItemCount
				% Load only timestamps and unit ID for selected channel
				[s(n), tmp.TimeStamp(m), ~, ~, tmp.UnitID(m)] = ...
						ns_GetSegmentData(Hfile, SegmentList(n), m);
			end			
		end
	end		
	d(n) = tmp; %#ok<AGROW>
	clear tmp;
end


