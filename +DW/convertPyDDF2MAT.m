function D = convertPyDDF2MAT(filename, varargin)
%------------------------------------------------------------------------
% D = convertPyDDF2MAT(varargin)
%------------------------------------------------------------------------
% loads text file output from PyDDF file in <filename>  
% stores data in struct and saves to a .mat file
%
% Each entity contains one or more indexed data entries that are ordered by
% increasing time.
% 
% The structure of the indexed data entries for each entity depends on the
% entity type:
% 
% List of EntityIDs needed to retrieve the information and data:
% 
% 	Event Entities (1)
% 		Discrete events that consist of small time-stamped text
% 		or binary data packets.  These are used to represent data such as trial
% 		markers, experimental events, digital input values, and embedded user
% 		comments.
% 		Each index of an event entity refers to a timestamp and data
% 		combination. The number of indexes is equal to the number of event
% 		entries for that event entity in the data file.
% 
% 	Analog Entities (2)
% 		Continuous, sampled data that represent digitized
% 		analog signals such as position, force, and other experiment signals, as
% 		well as electrode signals such as EKG, EEG and extracellular
% 		microelectrode recordings.
% 		Each index of an analog entity refers to a specific digitized sample.
% 		Each analog entity contains samples from a single channel and the
% 		number of indexes is equal to the number of samples present for that
% 		channel.
% 
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
% 	Neural Event Entities (4)
% 		Timestamps of event and segment entitities that
% 		are known to represent neural action potential firing times.  For
% 		example, if a segment entity contains sorted neural spike waveforms, each
% 		sorted unit is also exported as a neural entity.
% 		Each index of a neural event entity refers to a timestamp for each
% 		neural event.  Each neural event entity contains event times for a
% 		single neural source.  The number of indexes is equal to the number
% 		of entries for that neural event entity.
%------------------------------------------------------------------------
% Input Arguments:
%	filename		.ddf file from datawave (make sure path is included)
%
% Options (provided as '<option_name>', '<option_value>' pair):
%
%	'MATFILE', '<.output filename>'
%
% Options to load entity data (no option values needed):
% 	'EVENT'						loads event data
% 	'ANALOG'						loads analog data
% 	'SEGMENT'					loads segment data
% 	'NEURAL'						loads neural data
% 	'ALL'							loads all (event, analog segment, neural) data
%
% Output Arguments:
%	D			Data Structure
%------------------------------------------------------------------------
% See also: Neuroshare MATLAB API, DataWave Documentation, PyDDF
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 26 March, 2014 (SJS)
%
% Revisions:
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

%------------------------------------------------------------
%% read in file information header
%------------------------------------------------------------

if ~exist(filename, 'file')
	error('%s: input file %s not found', mfilename, filename);
end

fp = fopen(filename, 'r');

tag = fgetl(fp)

% check to make sure it's 'METADATA'
if ~strcmpi(tag, 'METADATA')
	fclose(fp);
	error('%s: no METADATA tag at header of file %s', mfilename, filename);
end

clear FileInfo
% read in metadata
for n = 1:14
	% read line
	tmpl = fgetl(fp);
	% parse into field:value
	[field, value] = readTaggedField(tmpl);
	FileInfo.(field) = value;
end

% convert entity count into number
FileInfo.EntityCount = str2num(FileInfo.EntityCount);

% read in entity info
tag = fgetl(fp)
% check to make sure it's 'ENTITY_INFO'
if ~strcmpi(tag, 'ENTITY_INFO')
	fclose(fp);
	error('%s: ENTITY_INFO not found after METADATA, file %s', mfilename, filename);
end
tmpl = fgetl(fp);
[tmp1, tmp2] = readTaggedField(tmpl);
Nentities = str2num(tmp2);
if Nentities ~= FileInfo.EntityCount
	fclose(fp);
	error('%s: Nentities ~= EntityCount, file %s', mfilename, filename);
end

% read entities
clear EntityInfo
for n = 1:FileInfo.EntityCount
	% read in line
	tmpl = fgetl(fp);
	% look for :
	semicols = find(tmpl == ':');
	% first three semicols are the relevant ones!
	% number/character after 1st semicolon is integer indicating
	% entity type
	EntityInfo(n).EntityType = str2num(tmpl(semicols(1)+1));
	% portion after second semicolon is entity label
	EntityInfo(n).EntityLabel = tmpl((semicols(3)+1):end);
	% store types
	entity_type_list(n) = EntityInfo(n).EntityType;
end

% build lists of types
EventList = find(entity_type_list == 1)
AnalogList = find(entity_type_list == 2)
SegmentList = find(entity_type_list == 3)
NeuralList = find(entity_type_list == 4)
% How many of a particular entity do we have
nNeural = length(NeuralList);       
nSegment = length(SegmentList);
nAnalog = length(AnalogList);
nEvent = length(EventList);
fprintf('%d Entities found\n', FileInfo.EntityCount);
fprintf('\t%d\tentity events available\n', nEvent);
fprintf('\t%d\tanalog events available\n', nAnalog);
fprintf('\t%d\tsegment events available\n', nSegment);
fprintf('\t%d\tneural events available\n', nNeural);

%--------------------------------------------------------------------------
% read R(ight) sound event markers
%--------------------------------------------------------------------------
tag = fgetl(fp)
% check to make sure it's 'MARKER'
if ~strcmpi(tag, 'MARKER')
	fclose(fp);
	error('%s: MARKER not found after ENTITY_INFO, file %s', mfilename, filename);
end
% get the label for this marker
tmpl = fgetl(fp);
[tmp1, tmp2] = readTaggedField(tmpl);
MarkerLabel = tmp2
% get the # of markers
tmpl = fgetl(fp);
[tmp1, tmp2] = readTaggedField(tmpl);
Rmarkers = str2num(tmp2)
% find corresponding EntityInfo member and store count
for n = 1:FileInfo.EntityCount
	if strcmp(EntityInfo(n).EntityLabel, MarkerLabel)
		EntityInfo(n).ItemCount = Rmarkers
	end
end
% get csv desc
tmpl = fgetl(fp);
[tmp1, CSVDesc] = readTaggedField(tmpl);

% pre allocate Event struct 
Event = repmat( struct('Info', [], 'EventCount', [], 'TimeStamp', [], ...
									'Data', {}, 'DataSize', []), 3, 1);
% this is a place holder
Event(1).Info.CSVDesc = 'null,null,null,';

% read in RIGHT markers, Store in Event list
Event(2).EventCount = Rmarkers;
Event(2).TimeStamp = zeros(1, Rmarkers);
Event(2).Data = cell(1, Rmarkers);
Event(2).DataSize = zeros(1, Rmarkers);
Event(2).Info.CSVDesc = CSVDesc;
for n = 1:Rmarkers
	tmpl = fgetl(fp);
	commaloc = find(tmpl == ',');
	% timestamp is first element
	Event(2).TimeStamp(n) = str2num(tmpl(1:(commaloc(1) - 1)));
	% Data{n} are from first comma to last
	Event(2).Data{n} = tmpl( (commaloc(1) + 1) : (commaloc(end) -1) );
end

%--------------------------------------------------------------------------
% read L(eft) sound event markers
%--------------------------------------------------------------------------
tag = fgetl(fp)
% check to make sure it's 'MARKER'
if ~strcmpi(tag, 'MARKER')
	fclose(fp);
	error('%s: MARKER not found after R MARKER, file %s', mfilename, filename);
end
% get the label for this marker
tmpl = fgetl(fp);
[tmp1, tmp2] = readTaggedField(tmpl);
MarkerLabel = tmp2
% get the # of markers
tmpl = fgetl(fp);
[tmp1, tmp2] = readTaggedField(tmpl);
Lmarkers = str2num(tmp2)
% find corresponding EntityInfo member and store count
for n = 1:FileInfo.EntityCount
	if strcmp(EntityInfo(n).EntityLabel, MarkerLabel)
		EntityInfo(n).ItemCount = Lmarkers
	end
end
% get csv desc
tmpl = fgetl(fp);
[tmp1, CSVDesc] = readTaggedField(tmpl);

% read in LEFT markers, Store in Event list
Event(3).EventCount = Lmarkers;
Event(3).TimeStamp = zeros(1, Lmarkers);
Event(3).Data = cell(1, Lmarkers);
Event(3).DataSize = zeros(1, Lmarkers);
Event(3).Info.CSVDesc = CSVDesc;
for n = 1:Rmarkers
	tmpl = fgetl(fp);
	commaloc = find(tmpl == ',');
	% timestamp is first element
	Event(3).TimeStamp(n) = str2num(tmpl(1:(commaloc(1) - 1)));
	% Data{n} are from first comma to last
	Event(3).Data{n} = tmpl( (commaloc(1) + 1) : (commaloc(end) -1) );
end


%--------------------------------------------------------------------------
% Segments
%--------------------------------------------------------------------------

% pre allocate Segment struct
Segment = repmat( struct('Info', [], 'SourceInfo', [], 'ItemCount', [], ...
									'TimeStamp', [], 'WaveForm', {}, ...
									'Nsamples', [], 'UnitID', []), nSegment, 1);

tag = fgetl(fp)
% check to make sure it's 'SEGMENT'
if ~strcmpi(tag, 'SEGMENT')
	fclose(fp);
	error('%s: SEGMENT not found after MARKER, file %s', mfilename, filename);
end
% get label and store it
tmpl = fgetl(fp);
[tmp1, label] = readTaggedField(tmpl);
Segment(1).SourceInfo.ProbeInfo = label;

% read max, min samples
[tmp1, tmp2] = readTaggedField(fgetl(fp));
Segment(1).Info.MaxSampleCount = str2num(tmp2);
[tmp1, tmp2] = readTaggedField(fgetl(fp));
Segment(1).Info.MinSampleCount = str2num(tmp2);
% read sample rate
[tmp1, tmp2] = readTaggedField(fgetl(fp));
Segment(1).Info.SampleRate = str2num(tmp2);
% read sample units
[tmp1, tmp2] = readTaggedField(fgetl(fp));
Segment(1).Info.Units = tmp2;

% store # of segments for this block
tmpl = fgetl(fp);
[tmp1, tmp2] = readTaggedField(tmpl);
Nsegments = str2num(tmp2)
% allocate things
Segment(1).ItemCount = Nsegments;
Segment(1).TimeStamp = zeros(Nsegments, 1);
Segment(1).WaveForm = cell(Nsegments, 1);
Segment(1).Nsamples = zeros(Nsegments, 1);
Segment(1).UnitID = zeros(Nsegments, 1);
% loop and read in data
for n = 1:Nsegments
	tmpl = fgetl(fp);
	T = sscanf(tmpl, '%f,');
	Segment(1).TimeStamp(n) = T(1);
	Segment(1).UnitID(n) = cast(T(2), 'int32');
	Segment(1).Nsamples(n) = cast(T(3), 'int32');
	Segment(1).WaveForm{n} = T(4:end)';
end

% check if we're done
if feof(fp) || (nSegment == 1)
	fclose(fp);
	if nSegment > 1
		warning('Only read 1 segment of data');
	end
	return
end

% continue with another segment
tag = fgetl(fp)
% check to make sure it's 'SEGMENT'
if ~strcmpi(tag, 'SEGMENT')
	fclose(fp);
	error('%s: SEGMENT not found after 1st SEGMENT, file %s', mfilename, filename);
end

% get label and store it
tmpl = fgetl(fp);
[tmp1, label] = readTaggedField(tmpl);
Segment(2).SourceInfo.ProbeInfo = label;
% read max, min samples
[tmp1, tmp2] = readTaggedField(fgetl(fp));
Segment(2).Info.MaxSampleCount = str2num(tmp2);
[tmp1, tmp2] = readTaggedField(fgetl(fp));
Segment(2).Info.MinSampleCount = str2num(tmp2);
% read sample rate
[tmp1, tmp2] = readTaggedField(fgetl(fp));
Segment(2).Info.SampleRate = str2num(tmp2);
% read sample units
[tmp1, tmp2] = readTaggedField(fgetl(fp));
Segment(2).Info.Units = tmp2;
% store # of segments for this block
tmpl = fgetl(fp);
[tmp1, tmp2] = readTaggedField(tmpl);
Nsegments = str2num(tmp2)
% allocate things
Segment(2).ItemCount = Nsegments;
Segment(2).TimeStamp = zeros(Nsegments, 1);
Segment(2).WaveForm = cell(Nsegments, 1);
Segment(2).Nsamples = zeros(Nsegments, 1);
Segment(2).UnitID = zeros(Nsegments, 1);
% loop and read in data
for n = 1:Nsegments
	tmpl = fgetl(fp);
	T = sscanf(tmpl, '%f,');
	Segment(2).TimeStamp(n) = T(1);
	Segment(2).UnitID(n) = cast(T(2), 'int32');
	Segment(2).Nsamples(n) = cast(T(3), 'int32');
	Segment(2).WaveForm{n} = T(4:end)';
end

fclose(fp);

%------------------------------------------------------------
% assign elements to D (output) struct
%------------------------------------------------------------
D = FileInfo;
D.EntityInfo = EntityInfo;
D.EventList = EventList;
D.AnalogList = AnalogList;
D.SegmentList = SegmentList;
D.NeuralList = NeuralList;
D.nNeural = nNeural;
D.nSegment = nSegment;
D.nAnalog = nAnalog;
D.nEvent = nEvent;

D.Event = Event;
D.Segment = Segment;

%------------------------------------------------------------
% clear unused data
%------------------------------------------------------------
clear FileInfo EntityInfo EventList AnalogList SegmentList NeuralList
clear nNeural nSegment nAnalog nEvent
clear Event
clear Segment



