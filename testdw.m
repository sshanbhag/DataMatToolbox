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
Events = D.loadMarkers;

%------------------------------------------------------------
%% parse markers
%------------------------------------------------------------
DataWaveDefaults

% create a list of Event IDS for each event entity.
%	EventID == 0 --> non-marker entity
%	EventID == 1 --> L channel markers
%	EventID == 2 --> R channel markers
EventID = zeros(D.DDF.nEvent, 1);

% first, find the R and L marker events (R is usually Events(2), 
% L is usually Events(3), but try not to assume)
for v = 1:D.DDF.nEvent
	% scan the CSVDesc from the NeuroShare Events.Info to get
	% a list (usually partial - for sound stimulus markers, the 
	% DataWave Neuroshare (or Matlab API, unknown at moment...) leaves
	% out some of the labels
	tmp = textscan(Events(v).Info.CSVDesc, '%s', ...
									'Delimiter', ',', ...
									'MultipleDelimsAsOne', 1);
	field_names = tmp{1};
	clear tmp;

	% check if the current Event entity has 'SoundType' as part of
	% the first element of its fields
	if strncmpi(field_names{1}, 'SoundType', length('SoundType'))
		% if so, assume that this is a "Marker" Entity and check for
		% R or L to determine Right or Left sound marker
		if strcmpi(field_names{1}(end), 'L')
			% left marker
			EventID(v) = L;
		elseif strcmpi(field_names{1}(end), 'R')
			% right marker
			EventID(v) = R;
		else
			EventID(v) = 0;
		end
	end
end

% check to ensure that both L and R Entities were found
if ~all( any(EventID == L) & any(EventID == R) )
	% if not abort
	error('%s: Markers for both (L and R) channels not found!', mfilename);
end

MarkerEvents = find(EventID);
nMarkerEvents = length(MarkerEvents);
if nMarkerEvents ~= 2
	error('%s: nMarkerEvents is %d (must be equal to 2)',...
					mfilename, nMarkerEvents);
end
% make sure # of events are the same
if ~(Events(MarkerEvents(1)).EventCount == Events(MarkerEvents(2)).EventCount)
	error('%s: EventCount mismatch!', mfilename);
else
	EventCount = Events(MarkerEvents(1)).EventCount; 
end

evL = find(EventID == L);
evR = find(EventID == R);
fprintf('Left Markers -> Event(%d)\n', evL);
fprintf('Right Markers -> Event(%d)\n', evR);

M = repmat(DW.Marker, EventCount, 1);

for n = 1:EventCount
	% get the right and left strings (Event Data)
	rstr = Events(evR).Data(n);
	lstr = Events(evL).Data(n);
	% parse into cell arrays, to preserve empty fields, do NOT
	% treat successive delimiters/whitespace as one (in call to csvscan)
	tmpR = csvscan(Events(evR).Data{n}, 0);
	tmpL = csvscan(Events(evL).Data{n}, 0);
	
	% ASSUME that one of the lists will not have an outputfile field.
	% pad with empty values 
	if length(tmpR) ~= MARKER_NBASE
		dlen = MARKER_NBASE - length(tmpR);
		if dlen < 0
			error('%d: length(tmpR) > MARKER_NBASE');
		else
			tmpR = [tmpR; cell(dlen, 1)]; %#ok<AGROW>
		end
	end
	if length(tmpL) ~= MARKER_NBASE
		dlen = MARKER_NBASE - length(tmpL);
		if dlen < 0
			error('%d: length(tmpL) > MARKER_NBASE');
		else
			tmpL = [tmpL; cell(dlen, 1)]; %#ok<AGROW>
		end
	end

	tmp = [tmpR; tmpL];
	elist = cell(size(tmp));
	% convert tags appropriately
	for t = 1:MARKER_NMARKERS
		% check for type
		switch MARKER_TYPES{t}
			case 'char'
				elist{t} = tmp{t};
			case 'float'
				elist{t} = str2double(tmp{t});
			case 'int'
				elist{t} = str2double(tmp{t});
			otherwise
				elist{t} = tmp{t};
		end
	end	% END t -> MARKERS_NMARKERS
	
	% assign val to object
	M(n).setValuesFromEventList(elist);
	clear tmpR tmpL tmp elist
end
