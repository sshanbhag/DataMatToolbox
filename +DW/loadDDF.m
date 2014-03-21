function [D, varargout] = loadDDF(filename, varargin)
%------------------------------------------------------------------------
% [D, H] = loadDDF(varargin)
%------------------------------------------------------------------------
% loads DDF file info using Neuroshare
%
% Note that SciWorks must be installed on the machine on which this
% function is run.  The SciWorks USB key is NOT needed.
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
%	'DLL',	'<dll_name>'	.dll file to use for Neuroshare interface
%									default value is 'C:\DataWave\DWShared\nsDWFile.dll'
% 	'EVENT'						loads event data
% 	'ANALOG'						loads analog data
% 	'SEGMENT'					loads segment data
% 	'NEURAL'						loads neural data
% 	'ALL'							loads all (event, analog segment, neural) data
%
% Output Arguments:
%	D			Data Structure
%	H			handle to Neuroshare interface
%------------------------------------------------------------------------
% See also: Neuroshare MATLAB API
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 19 December, 2011 (SJS)
%
% Revisions:
%	3 Jan 2013 (SJS):
% 	 -	updated docs
% 	 -	added input options for loading data
% 	 - reworked outputs
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% some definitions
%------------------------------------------------------------------------
%------------------------------------------------------------------------
H = 0;
% default DLL path and name
DLLName = 'C:\DataWave\DWShared\nsDWFile.dll';
% Flags for loading data
EVENT = 0;
ANALOG = 0;
SEGMENT = 0;
NEURAL = 0;

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% check output arguments
%------------------------------------------------------------------------
%------------------------------------------------------------------------
if nargout > 1
	varargout{1} = [];
end
D = [];

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% check input arguments
%------------------------------------------------------------------------
%------------------------------------------------------------------------
if nargin == 0
	error('%s: need input filename\n', mfilename);
end
nvararg = length(varargin);
if nvararg
	aindex = 1;
	while aindex <= nvararg
		switch(upper(varargin{aindex}))
			% select DLL file
			case 'DLL'
				DLLName = varargin{aindex + 1};
				aindex = aindex + 2;
			% Load all data
			case 'ALL'
				EVENT = 1;
				ANALOG = 1;
				SEGMENT = 1;
				NEURAL = 1;
				aindex = aindex + 1;
			% load event
			case 'EVENT'
				EVENT = 1;
				aindex = aindex + 1;
			% Load analog data
			case 'ANALOG'
				ANALOG = 1;
				aindex = aindex + 1;
			% Load SEGMENT data
			case 'SEGMENT'
				SEGMENT = 1;
				aindex = aindex + 1;
			% Load NEURAL data
			case 'NEURAL'
				NEURAL = 1;
				aindex = aindex + 1;		
			% trap unknown options
			otherwise
				error('%s: Unknown option %s', mfilename, varargin{aindex});
		end		% END SWITCH
	end		% END WHILE aindex
end		% END IF nvararg
%------------------------------------------------------------------------
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% Load the appropriate DLL
%------------------------------------------------------------------------
%------------------------------------------------------------------------
[nsresult] = ns_SetLibrary(DLLName);
if (nsresult ~= 0)
	fprintf('%s: DLL %s was not found or load error!\n\n', mfilename, DLLName);
	return
else
	fprintf('\t%s: DLL loaded\n', mfilename);
	% get information about library
	[ns_result, nsLibraryInfo] = ns_GetLibraryInfo();
	disp(nsLibraryInfo)
end
%------------------------------------------------------------------------
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% Load data file and display some info about the file
%------------------------------------------------------------------------
%------------------------------------------------------------------------
[nsresult, H] = ns_OpenFile(filename);
if (nsresult ~= 0)
	fprintf('%s: Data file %s did not open!\n', mfilename, filename);
	return
else
	% Get file information
	[nsresult, FileInfo] = ns_GetFileInfo(H);
	% Gives you EntityCount, TimeStampResolution and TimeSpan
	if (nsresult ~= 0)
		fprintf('%s: Data file %s information did not load!\n', mfilename, filename);
		return
	else
		fprintf('\n\nFileInfo:\n');
		disp(FileInfo);
	end
end
%------------------------------------------------------------------------
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% Build catalogue of entities
%------------------------------------------------------------------------
%------------------------------------------------------------------------
if FileInfo.EntityCount
	[nsresult, EntityInfo] = ns_GetEntityInfo(H, (1:FileInfo.EntityCount));
else
	% Close data file. Should be done by the library but just in case. 
	ns_CloseFile(H);
	error('%s: no Entities found in file %s', mfilename, filename);
end
EventList = find([EntityInfo.EntityType] == 1);
AnalogList = find([EntityInfo.EntityType] == 2);
SegmentList = find([EntityInfo.EntityType] == 3);
NeuralList = find([EntityInfo.EntityType] == 4);

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

%------------------------------------------------------------
% clear unused data
%------------------------------------------------------------
clear FileInfo EntityInfo EventList AnalogList SegmentList NeuralList
clear nNeural nSegment nAnalog nEvent

%------------------------------------------------------------
%% About Entities
%------------------------------------------------------------
% Each entity contains one or more indexed data entries that are ordered by
% increasing time.  The API provides functions for querying the
% characteristics of the file, the number of entities, and the
% characteristics of each entity, including the number of indexes for each
% entity.
% 
% The structure of the indexed data entries for each entity depends on the
% entity type (Event, Analog, Segment, Neural) and, as a
% result, they will be read from the file in different ways
%------------------------------------------------------------

%------------------------------------------------------------
%% Read Event Entities
%------------------------------------------------------------
% Each index of an event entity refers to a timestamp and data combination.
% The number of indexes is equal to the number of event entries for that
% event entity in the data file.
%------------------------------------------------------------
if EVENT
	D.Event = readEvent(H, D);
end
%------------------------------------------------------------
%% read analog entities
%------------------------------------------------------------
% Each index of an analog entity refers to a specific digitized sample.
% Each analog entity contains samples from a single channel and the number
% of indexes is equal to the number of samples present for that channel.
%------------------------------------------------------------
if ANALOG
	D.Analog = getAnalog(H, D);
end
%------------------------------------------------------------
%% read segment entities
%------------------------------------------------------------
% Each index of a segment entity refers to a short, time-stamped segment of
% analog data from one or more sources.  The number of indexes is equal to
% the number of entries for that segment entity in the file.
%------------------------------------------------------------
if SEGMENT
	D.Segment = readSegment(H, D);
end
%------------------------------------------------------------
%% read neural entitites (if they exist)
%------------------------------------------------------------
% Each index of a neural event entity refers to a timestamp for each neural
% event.  Each neural event entity contains event times for a single neural
% source.  The number of indexes is equal to the number of entries for that
% neural event entity.
%------------------------------------------------------------
if NEURAL
	D.Neural = readNeural(H, D);
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% build output args
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% if output for file handle object was requested, assign it here, 
if nargout > 1
	varargout{1} = H;
% otherwise, close the file
else
	ns_CloseFile(H);
end
