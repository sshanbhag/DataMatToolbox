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

for n = 1:EventCount
	% get the right and left strings (Event Data)
	rstr = Events(evR).Data(n);
	lstr = Events(evL).Data(n);
	% parse into cell arrays
	tmpR = csvscan(Events(evR).Data{n}, 0);
	tmpL = csvscan(Events(evL).Data{n}, 0);
keyboard
	M = DW.Marker('EVENT_STRING', [tmpR{1}; tmpL{1}]);

	clear tmpR tmpL
end

field_names = cell(D.DDF.nEvent, 1);
field_vals = cell(D.DDF.nEvent, 1);
for v = 1:D.DDF.nEvent
	tmp = textscan(Events(v).Info.CSVDesc', '%s', ...
									'Delimiter', ',', ...
									'MultipleDelimsAsOne', 1);
	field_names{v} = tmp{1};
	field_vals{v} = cell(Events(v).EventCount, 1);	
	if Events(v).EventCount
		for n = 1:Events(v).EventCount
			tmp = textscan(Events(v).Data{n}', '%s', ...
										'Delimiter', ',', ...
										'MultipleDelimsAsOne', 1);
			field_vals{v}{n} = tmp{1};
		end	% END n
		

	end	% END if EventCount
end	% END v DDF.nEvent index

for v = 1:length(field_names)
	fprintf('event %d:\n', v)
	disp(field_names{v})
	if ~isempty(field_vals{v})
		disp(field_vals{v}{1})
	end
	fprintf('\n\n');
end
