%-----------------------------------------------------------------------------
% NS.m
%-----------------------------------------------------------------------------
% DataMat Toolbox
% DW Package
% Class Definition
%-----------------------------------------------------------------------------
%	NS facilitates access to NeuroShare data
%
%  obj = DWdata('xxx_y_zz_ttt_n.mat') loads the following:
%     obj.fpath;   %path to file
%     obj.fname;   %name of file
%     obj.fext;    %file extension
%     obj.Background; %really this should be a class too...
%     obj.D;          %really this should be a class too...
%     obj.Stimulus;   %really this should be a class too...
%     obj.fullfname;  %full file name
%-----------------------------------------------------------------------------
% See also: NeuroShare Toolbox for Matlab
%-----------------------------------------------------------------------------

%-----------------------------------------------------------------------------
%	Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%-----------------------------------------------------------------------------
% Created: 9 January, 2013 (SJS)
%
% Revisions:
%	25 Jan 2013 (SJS): added some methods to get segment info
%-----------------------------------------------------------------------------
% TO DO: lots....
%-----------------------------------------------------------------------------

%*****************************************************************************
%*****************************************************************************
%*****************************************************************************
% class definition
%*****************************************************************************
%*****************************************************************************
%*****************************************************************************
classdef (ConstructOnLoad = true) NS < handle
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Define protected properties
	%------------------------------------------------------------------------
	properties (SetAccess = protected)
		% Library
		nsLibraryInfo;
		libraryLoaded = 0;
		
		% data and information
		FileType;
		EntityCount;
		TimeStampResolution;
		TimeSpan;
		AppName;
		Time_Year;
		Time_Month;
		Time_Day;
		Time_Hour;
		Time_Min;
		Time_Sec;
		Time_MilliSec;
		FileComment;
		EntityInfo;
		EventList;
		AnalogList;
		SegmentList;
		NeuralList;
		nNeural;
		nSegment;
		nAnalog;
		nEvent;
		% neuroshare file handle
		Hfile;
		% data file
		datapath;
		dataname;
		dataext;
		% neuroshare DLL information
		DLLpath;
		DLLname;
		DLLfullname;
	end
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------

	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Define other properties
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	properties
		datafullname;  %full file name
	end
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
  
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Define PUBLIC methods
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	methods
	
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		% Constructor
		%------------------------------------------------------------------------
		function obj = NS(varargin)
		%------------------------------------------------------------------------
		% obj = NS(filename, varargin)
		%------------------------------------------------------------------------
		% Initializes the object
		%------------------------------------------------------------------------
		% Input Arguments:
		%	filename		.ddf file from datawave (make sure path is included)
		%
		% Options (provided as '<option_name>', '<option_value>' pair):
		%	'DLL',	'<dll_name>'	.dll file to use for Neuroshare interface
		%										default value is 
		%										'C:\DataWave\DWShared\nsDWFile.dll'
		% 	'EVENT'						loads event data
		% 	'ANALOG'						loads analog data
		% 	'SEGMENT'					loads segment data
		% 	'NEURAL'						loads neural data
		% 	'ALL'							loads all (event, analog segment, neural) data
		%
		% Output Arguments:
		%	obj		NS object
		%------------------------------------------------------------------------
		% See also: Neuroshare MATLAB API
		%------------------------------------------------------------------------

			%--------------------
			% set defaults
			%--------------------
			% default DLL path and name
			obj.DLLpath = 'C:\DataWave\DWShared';
			obj.DLLname = 'nsDWFile.dll';
			obj.DLLfullname = fullfile(obj.DLLpath, obj.DLLname);
			% Flags for loading data
			EVENT = 0;
			ANALOG = 0;
			SEGMENT = 0;
			NEURAL = 0;
			%--------------------
			% parse input args
			%--------------------
			% get filename information - this must be first input arg
			% (required)!!
			if nargin == 0
				error('%s: need input filename\n', mfilename);
			else
				[obj.datapath, obj.dataname, obj.dataext] = fileparts(varargin{1});
				obj.datafullname = ...
									fullfile(obj.datapath, [obj.dataname obj.dataext]);
				if ~exist(obj.datafullname, 'file')
					error('%s: datafile %s not found!', mfilename, obj.datafullname);
				end
			end
			arg = 2;
			while arg <= nargin
				switch(upper(varargin{arg}))
					% select DLL file
					case 'DLL'
						[obj.DLLpath, obj.DLLname] = varargin{arg + 1};
						obj.DLLfullname = fullfile(obj.DLLpath, obj.DLLname);
						if ~exist(obj.DLLfullname, 'file')
							error('%s: DLL %s not found!', mfilename, obj.DLLfullname);
						end
						arg = arg + 2;
					% Load all data
					case 'ALL'
						EVENT = 1;
						ANALOG = 1;
						SEGMENT = 1;
						NEURAL = 1;
						arg = arg + 1;
					% load event
					case 'EVENT'
						EVENT = 1;
						arg = arg + 1;
					% Load analog data
					case 'ANALOG'
						ANALOG = 1;
						arg = arg + 1;
					% Load SEGMENT data
					case 'SEGMENT'
						SEGMENT = 1;
						arg = arg + 1;
					% Load NEURAL data
					case 'NEURAL'
						NEURAL = 1;
						arg = arg + 1;
					% trap unknown options
					otherwise
						error('%s: Unknown option %s', mfilename, varargin{arg});
				end		% END SWITCH
			end		% END WHILE arg
			
			%-----------------------------
			% set the NeuroShare Library
			%-----------------------------
			obj.setLibrary;
			%-----------------------------
			% get file information
			%-----------------------------
			obj.getFileInfo;
		end	% end NS (constructor)
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function obj = setLibrary(obj)
		%------------------------------------------------------------------------
		% NS.setLibrary
		%------------------------------------------------------------------------
		% sets NeuroShare DLL and loads library information
		%------------------------------------------------------------------------
			nsresult = ns_SetLibrary(obj.DLLfullname);
			if (nsresult ~= 0)
				fprintf('%s: DLL %s load error!\n\n', mfilename, obj.DLLfullname);
				return
			else
				fprintf('%s:\tDLL %s loaded\n', mfilename, obj.DLLfullname);
				% get information about library
				fprintf('\tnsLibraryInfo:\n');
				[ns_result, obj.nsLibraryInfo] = ns_GetLibraryInfo();
				obj.libraryLoaded = 1;
				disp(obj.nsLibraryInfo)
			end
		end	% END setLibrary
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------

		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function obj = getFileInfo(obj)
		%------------------------------------------------------------------------
		% NS.getFileInfo(obj)
		%------------------------------------------------------------------------
		% get some info about the file
		%------------------------------------------------------------------------
		
			%----------------------------------------------------
			% open file using NeuroShare interface
			%----------------------------------------------------
			obj.openFile;
			%----------------------------------------------------
			% Get file information: EntityCount, TimeStampResolution and TimeSpan
			%----------------------------------------------------
			[nsresult, FileInfo] = ns_GetFileInfo(obj.Hfile);
			if (nsresult ~= 0)
				fprintf('%s: file information did not load!\n', mfilename);
				obj.closeFile;
				return
			end
			disp(FileInfo);
			%----------------------------------------------------
			% assign field data from FileInfo to corresponding object params
			%----------------------------------------------------
			fnames = fieldnames(FileInfo);
			for n = 1:length(fnames)
				obj.(fnames{n}) = FileInfo.(fnames{n});
			end
			%----------------------------------------------------
			% Build catalogue of entities if any exist
			%----------------------------------------------------
			if FileInfo.EntityCount
				[nsresult, obj.EntityInfo] = ...
								ns_GetEntityInfo(obj.Hfile, (1:FileInfo.EntityCount));
				for n = 1:FileInfo.EntityCount
					switch(obj.EntityInfo(n).EntityType)
						case 1
							obj.EntityInfo(n).EntityTypeStr = 'EVENT';
						case 2
							obj.EntityInfo(n).EntityTypeStr = 'ANALOG';
						case 3
							obj.EntityInfo(n).EntityTypeStr = 'SEGMENT';
						case 4
							obj.EntityInfo(n).EntityTypeStr = 'NEURAL';
						otherwise
							obj.EntityInfo(n).EntityTypeStr = 'UNKNOWN';
					end	
				end
			else
				fprintf('%s: no Entities found in file %s', ...
								mfilename, obj.datafullname);
				return
			end
			obj.EventList = find([obj.EntityInfo.EntityType] == 1);
			obj.AnalogList = find([obj.EntityInfo.EntityType] == 2);
			obj.SegmentList = find([obj.EntityInfo.EntityType] == 3);
			obj.NeuralList = find([obj.EntityInfo.EntityType] == 4);
			%----------------------------------------------------
			% How many of a particular entity do we have
			%----------------------------------------------------
			obj.nNeural = length(obj.NeuralList);       
			obj.nSegment = length(obj.SegmentList);
			obj.nAnalog = length(obj.AnalogList);
			obj.nEvent = length(obj.EventList);
			%----------------------------------------------------
			% give info to user
			%----------------------------------------------------
			fprintf('%d Entities found\n', FileInfo.EntityCount);
			fprintf('\t%d\tentity events available\n', obj.nEvent);
			fprintf('\t%d\tanalog events available\n', obj.nAnalog);
			fprintf('\t%d\tsegment events available\n', obj.nSegment);
			fprintf('\t%d\tneural events available\n', obj.nNeural);
			%----------------------------------------------------
			% clear unused data
			%----------------------------------------------------
			clear FileInfo
			%----------------------------------------------------			
			% close file
			%----------------------------------------------------
			obj.closeFile;
		end	% END getFileInfo
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------

		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		% About Entities
		%------------------------------------------------------------
		% Each entity contains one or more indexed data entries that are ordered
		% by increasing time.  The API provides functions for querying the
		% characteristics of the file, the number of entities, and the
		% characteristics of each entity, including the number of indexes for each
		% entity.
		% 
		% The structure of the indexed data entries for each entity depends on the
		% entity type (Event, Analog, Segment, Neural) and, as a
		% result, they will be read from the file in different ways
		%------------------------------------------------------------		
		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function Events = getEvents(obj)
		%------------------------------------------------------------------------
		% Events = NS.getEvents
		%------------------------------------------------------------------------
		% Read Event Entities
		%------------------------------------------------------------------------
		% 	Event Entities (entity type = 1)
		% 	Discrete events that consist of small time-stamped text
		% 	or binary data packets.  These are used to represent data such as trial
		% 	markers, experimental events, digital input values, and embedded user
		% 	comments.
		% 	Each index of an event entity refers to a timestamp and data
		% 	combination. The number of indexes is equal to the number of event
		% 	entries for that event entity in the data file.
		%------------------------------------------------------------------------

			%----------------------------------------------------
			% setup
			%----------------------------------------------------
			Events = [];
			if (obj.nEvent == 0)
				disp('No event entities available!');
				return
			end
			obj.openFile;
			
			%----------------------------------------------------
			% preallocate events
			%----------------------------------------------------
			Events = repmat( struct(	'Info', [], ...
												'EventCount', [], ...
												'TimeStamp', [], ...
												'Data', {}, ...
												'DataSize', []	), ...
									obj.nEvent, 1	);
			s = zeros(obj.nEvent, 1);
			%----------------------------------------------------
			% Read Events in the eventlist
			%----------------------------------------------------
			for n = 1:obj.nEvent
				% get event information
				[s(n, 1), Events(n).Info] = ...
											ns_GetEventInfo(obj.Hfile, obj.EventList(n));
				% # of Items for this event
				 Events(n).EventCount = obj.EntityInfo(obj.EventList(n)).ItemCount;
				if Events(n).EventCount
					% read data for this event if length ~= 0
					% pre-allocate arrays
					s2 = zeros(Events(n).EventCount, 1);
					Events(n).TimeStamp = zeros(size(s2));
					Events(n).Data = cell(size(s2));
					Events(n).DataSize = zeros(size(s2));
					% read events
					for m = 1:Events(n).EventCount
						[s2(m), ...
						 Events(n).TimeStamp(m), ...
						 tmp, ...
						 Events(n).DataSize(m)] = ...
										ns_GetEventData(obj.Hfile, obj.EventList(n), m);
						% ns_GetEventData() returns a cell array of size {1} for Data.
						% store the contents in the {m, 1} cell array in Events.
						Events(n).Data{m} = tmp{1};
					end
					clear tmp s2
				else
					% otherwise, assign empty vals to data structure
					Events(n).TimeStamp = [];
					Events(n).Data = {};
					Events(n).DataSize = [];
				end
			end	% END nEvent loop
			%----------------------------------------------------
			% close file
			%----------------------------------------------------
			obj.closeFile;
		end	% END getEvents
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------

		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function Analog = getAnalog(obj, varargin)
		%------------------------------------------------------------------------
		% Analog = NS.getAnalog(varargin)
		%------------------------------------------------------------------------
		% Read analog entities
		%------------------------------------------------------------------------
		% 	Analog Entities (2)
		% 	Continuous, sampled data that represent digitized
		% 	analog signals such as position, force, and other experiment signals,
		% 	as well as electrode signals such as EKG, EEG and extracellular
		% 	microelectrode recordings.
		% 	Each index of an analog entity refers to a specific digitized sample.
		% 	Each analog entity contains samples from a single channel and the
		% 	number of indexes is equal to the number of samples present for that
		% 	channel.
		%------------------------------------------------------------------------
		% There are two ways that are used for reading analog data
		% Analog = obj.getAnalog(<entity list (optional)>)
		% 	if no entity list (values from N.AnalogList) is provided,
		% 	all analog entities will be read.
		% 	Warning: This might use up all memory if the recordings are
		% 	of a long duration!
		% 
		% Analog = obj.getAnalog(<entity list>, <start point>, < points to read>)
		%------------------------------------------------------------------------
		
			%----------------------------------------------------
			% checks & setup
			%----------------------------------------------------
			if (obj.nAnalog == 0)
				disp('No analog entities available!');
				return
			end
			Analog = [];
			obj.openFile;
			%----------------------------------------------------
			% check if user asked for a specific list of analog entities
			% and/or a specific chunk of the waveform
			% note that even if length(varargin) is 0, nargin will be 1
			%----------------------------------------------------
			if nargin == 1
				anaList = [];
				anaPoints = {};
			elseif nargin == 2
				anaList = varargin{1};
			elseif nargin == 4
				anaList = varargin{1};
				anaPoints = {varargin{2}, varargin{3}};
			else
				fprintf('NS.getAnalog error: incorrect input args\n');
				Analog = [];
				return
			end
			% if anaList is empty, use the full AnalogList
			if isempty(anaList)
				anaList = obj.AnalogList;
			end
			%----------------------------------------------------
			% pre-allocate Analog and s arrays
			%----------------------------------------------------
			Nlist = length(anaList);
			Analog = repmat(	struct(	'Info',		[], ...
												'startPoint',	[], ...
												'nPoints',	[], ...
												'Data',		[] ), ...
									Nlist, 1	);
			s = zeros(Nlist, 1);
			% if anaPoints is empty
			if isempty(anaPoints)
				% use all points
				anaPoints{1} = ones(Nlist, 1);
				anaPoints{2} = zeros(Nlist, 1);
				for n = 1:Nlist
					anaPoints{2}(n) = obj.EntityInfo(anaList(n)).ItemCount;
				end
			end
			%----------------------------------------------------
			% read analog entity 
			%----------------------------------------------------
			for n = 1:Nlist
				% get the information about the current analog entity
				[s(n), Analog(n).Info] = ns_GetAnalogInfo(obj.Hfile, anaList(n));
				% read analog entity data
				[s(n), nRead, Analog(n).Data] = ...
								ns_GetAnalogData(	obj.Hfile, ...
														anaList(n), ...
														anaPoints{1}(n), ...
														anaPoints{2}(n)	);
				% store # of analog points (items) for this entity
				Analog(n).nPoints = length(Analog(n).Data);
				if Analog(n).nPoints ~= anaPoints{2}(n)
					fprintf('%s warn: nPoints ~= points requested\n', mfilename);
				end
				Analog(n).startPoint = anaPoints{1}(n);	
			end
			%----------------------------------------------------
			% close file
			%----------------------------------------------------
			obj.closeFile
		end	% END getAnalog
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------

		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function S = getSegment(obj, varargin)
		%------------------------------------------------------------------------
		% Segment = NS.getSegment(segList, wavList)
		%------------------------------------------------------------------------
		% read segment entities
		%------------------------------------------------------------------------
		% Segment Entities (3)
		% 	Short, time-stamped segments of digitized analog
		% 	signals in which the segments are separated by variable amounts of time.
		% 	Segment Entities can contain data from more than one source. They are
		% 	intended to represent discontinuous analog signals such as extracellular
		% 	spike waveforms from electrodes or groups of electrodes.
		% 	Each index of a segment entity refers to a short, time-stamped
		% 	segment of analog data from one or more sources.  The number of
		% 	indexes is equal to the number of entries for that segment entity in
		% 	the file.
		%------------------------------------------------------------------------
		%	segList: a list of segments to get (from NS.SegmentList)
		%				if not provided, all segments will be retrieved
		%					
		%	wavList:	a cell vector listing the waveforms to retrieve
		%				if not provided all waveforms will be retrieved and 
		%				you will be sad.
		%
		%				cell must be same length as segList (or NS.SegmentList)
		%				and will have format like:
		%				{ [indices of waves for SegmentList(1)], 
		%				  [indices of waves for SegmentList(2)], 
		%					  .
		%				     .
		%					[indices of waves for SegmentList(n)] }
		%
		%				To indicate that waveforms for a particular 
		%				segment will not be read, enter empty matrix, []
		% 
		% 	A zero unit ID is unclassified, then follow unit 1, 2, 3, etc. Unit
		% 	255 is noise.
		% 
		%------------------------------------------------------------------------

			%----------------------------------------------------
			% checks & setup
			%----------------------------------------------------
			if (obj.nSegment == 0)
				disp('No segment entities available!');
				return
			end
			S = [];
			
			%----------------------------------------------------
			% check if user asked for a specific list of analog entities
			% note that even if length(varargin) is 0, nargin will be 1
			%----------------------------------------------------
			if nargin == 1
				segList = [];
				wavList = {};
			elseif nargin == 2
				segList = varargin{1};
			elseif nargin == 3
				segList = varargin{1};
				wavList = varargin{2};
			else
				fprintf('NS.getSegment error: incorrect input args\n');
				S = [];
				return
			end
			%----------------------------------------------------
			% open file
			%----------------------------------------------------
			obj.openFile;
			
			%----------------------------------------------------
			% if segList is empty, use the full SegmentList
			%----------------------------------------------------
			if isempty(segList)
				segList = obj.SegmentList;
			end
			Nseg = length(segList);
			%----------------------------------------------------
			% if wavList is empty, use all waveforms for SegmentList
			% indicate this by creating cell vector of empty arrays
			%----------------------------------------------------
			if isempty(wavList)
				wavList = cell(Nseg, 1);
				for n = 1:Nseg
					wavList{n} = [];
				end
			end

			%----------------------------------------------------
			% pre-allocate things
			%----------------------------------------------------
			t = zeros(Nseg, 1);
			S = repmat(	struct(	'Info',			[], ...
										'SourceInfo',	[], ...
										'MaxCount',		[], ...
										'ItemCount',	[], ...
										'TimeStamp',	[], ...
										'WaveForm',		{}, ...
										'Nsamples',		[], ...
										'UnitID',		[]), ...
							Nseg, 1);

			%----------------------------------------------------
			% loop through segments
			%----------------------------------------------------
			for n = 1:Nseg
				% get info and source info
				[t(n), S(n).Info] = ...
										ns_GetSegmentInfo(obj.Hfile, segList(n));
				[t(n), S(n).SourceInfo] = ...
							ns_GetSegmentSourceInfo(obj.Hfile, segList(n), 1);
				% get # of items (waveforms) available for this segment
				S(n).MaxCount = obj.EntityInfo(segList(n)).ItemCount;
				% only load segments if MaxCount > 0
				if isempty(S(n).MaxCount) || (S(n).MaxCount == 0)
					fprintf('%s warning: no segment items available for %d\n', ...
														mfilename, segList(n))
				else
					if isempty(wavList{n})
						fprintf('%s: reading all %d waveforms for Segment %d\n', ...
										mfilename, S(n).MaxCount, n);
						wavList{n} = 1:S(n).MaxCount; %#ok<AGROW>
					end
					% get # of waves to read in.
					Nwav = length(wavList{n});
					% store number of wavs to be read in
					S(n).ItemCount = Nwav;
					% if # waves ~= 0, read them in
					if Nwav > 0
						% pre-allocate storage
						S(n).TimeStamp = zeros(Nwav, 1);
						S(n).WaveForm = cell(Nwav, 1);
						S(n).Nsamples = zeros(Nwav, 1);
						S(n).UnitID = zeros(Nwav, 1);
						for m = 1:Nwav
							% Load the waveforms on each selected channel
							[t(n), S(n).TimeStamp(m), ...
									 S(n).WaveForm{m}, ...
									 S(n).Nsamples(m), ...
									 S(n).UnitID(m) ] = ...
														ns_GetSegmentData(obj.Hfile, ...
																				segList(n), ...
																				wavList{n}(m) );
						end	% END m Nwav
					end	% END if Nwav
				end	% END if
			end	% END n
			%----------------------------------------------------
			% close file
			%----------------------------------------------------
			obj.closeFile;
		end	% END getSegment
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------

		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function S = getSegmentTimeStampsOnly(obj, varargin)
		%------------------------------------------------------------------------
		% Segment = NS.getSegmentTimeStampsOnly(segList, wavList)
		%------------------------------------------------------------------------
		% read segment entities
		%------------------------------------------------------------------------
		% Segment Entities (3)
		% 	Short, time-stamped segments of digitized analog
		% 	signals in which the segments are separated by variable amounts of time.
		% 	Segment Entities can contain data from more than one source. They are
		% 	intended to represent discontinuous analog signals such as extracellular
		% 	spike waveforms from electrodes or groups of electrodes.
		% 	Each index of a segment entity refers to a short, time-stamped
		% 	segment of analog data from one or more sources.  The number of
		% 	indexes is equal to the number of entries for that segment entity in
		% 	the file.
		%------------------------------------------------------------------------
		%	segList: a list of segments to get (from NS.SegmentList)
		%				if not provided, all segments will be retrieved
		%					
		%	wavList:	a cell vector listing the waveforms to retrieve
		%				if not provided all waveforms will be retrieved and 
		%				you will be sad.
		%
		%				cell must be same length as segList (or NS.SegmentList)
		%				and will have format like:
		%				{ [indices of waves for SegmentList(1)], 
		%				  [indices of waves for SegmentList(2)], 
		%					  .
		%				     .
		%					[indices of waves for SegmentList(n)] }
		%
		%				To indicate that waveforms for a particular 
		%				segment will not be read, enter empty matrix, []
		% 
		% 	A zero unit ID is unclassified, then follow unit 1, 2, 3, etc. Unit
		% 	255 is noise.
		% 
		%------------------------------------------------------------------------

			%----------------------------------------------------
			% checks & setup
			%----------------------------------------------------
			if (obj.nSegment == 0)
				disp('No segment entities available!');
				return
			end
			S = [];
			
			%----------------------------------------------------
			% check if user asked for a specific list of analog entities
			% note that even if length(varargin) is 0, nargin will be 1
			%----------------------------------------------------
			if nargin == 1
				segList = [];
				wavList = {};
			elseif nargin == 2
				segList = varargin{1};
			elseif nargin == 3
				segList = varargin{1};
				wavList = varargin{2};
			else
				fprintf('NS.getSegment error: incorrect input args\n');
				S = [];
				return
			end
			%----------------------------------------------------
			% open file
			%----------------------------------------------------
			obj.openFile;
			
			%----------------------------------------------------
			% if segList is empty, use the full SegmentList
			%----------------------------------------------------
			if isempty(segList)
				segList = obj.SegmentList;
			end
			Nseg = length(segList);
			%----------------------------------------------------
			% if wavList is empty, use all waveforms for SegmentList
			% indicate this by creating cell vector of empty arrays
			%----------------------------------------------------
			if isempty(wavList)
				wavList = cell(Nseg, 1);
				for n = 1:Nseg
					wavList{n} = [];
				end
			end

			% pre-allocate things
			t = zeros(Nseg, 1);
			S = repmat(	struct(	'Info',			[], ...
										'SourceInfo',	[], ...
										'MaxCount',		[], ...
										'ItemCount',	[], ...
										'TimeStamp',	[], ...
										'UnitID',		[]), ...
							Nseg, 1);

			% loop through segments
			for n = 1:Nseg
				% get info and source info
				[t(n), S(n).Info] = ...
										ns_GetSegmentInfo(obj.Hfile, segList(n));
				[t(n), S(n).SourceInfo] = ...
							ns_GetSegmentSourceInfo(obj.Hfile, segList(n), 1);
				% get # of items (waveforms) available for this segment
				S(n).MaxCount = obj.EntityInfo(segList(n)).ItemCount;
				% only load segments if MaxCount > 0
				if isempty(S(n).MaxCount) || (S(n).MaxCount == 0)
					fprintf('%s warning: no segment items available for %d\n', ...
														mfilename, segList(n))
				else
					if isempty(wavList{n})
						fprintf('%s: reading all %d waveforms for Segment %d\n', ...
										mfilename, S(n).MaxCount, n);
						wavList{n} = 1:S(n).MaxCount; %#ok<AGROW>
					end
					% get # of waves to read in.
					Nwav = length(wavList{n});
					% store number of wavs to be read in
					S(n).ItemCount = Nwav;
					% if # waves ~= 0, read them in
					if Nwav > 0
						% pre-allocate storage
						S(n).TimeStamp = zeros(Nwav, 1);
						S(n).WaveForm = cell(Nwav, 1);
						S(n).Nsamples = zeros(Nwav, 1);
						S(n).UnitID = zeros(Nwav, 1);
						for m = 1:Nwav
							% Load the timestamps on each selected channel
							[~, S(n).TimeStamp(m), ...
									 ~, ...
									 ~, ...
									 S(n).UnitID(m) ] = ...
														ns_GetSegmentData(obj.Hfile, ...
																				segList(n), ...
																				wavList{n}(m) );
	
						end	% END m
					end	% END if Nwav
				end	% END if
			end	% END n
			%----------------------------------------------------
			% close file
			%----------------------------------------------------
			obj.closeFile;
			clear t;
		end	% END getSegment
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function [SegmentInfo, SourceInfo] = getSegmentInfo(obj, varargin)
		%------------------------------------------------------------------------
		% [SegmentInfo, SourceInfo] = NS.getSegmentInfo(segList)
		%------------------------------------------------------------------------
		% read segment information
		%------------------------------------------------------------------------
		%	segList: a list of segments to get (from NS.SegmentList)
		%				if not provided, all segments will be retrieved
		%					
		%------------------------------------------------------------------------

			%----------------------------------------------------
			% checks
			%----------------------------------------------------
			if (obj.nSegment == 0)
				SegmentInfo = [];
				SourceInfo = [];
				disp('No segment entities available!');
				return
			end
			%----------------------------------------------------
			% check if user asked for a specific list of segment entities
			% note that even if length(varargin) is 0, nargin will be 1
			%----------------------------------------------------
			if nargin == 1
				segList = [];
			elseif nargin == 2
				segList = varargin{1};
			else
				fprintf('NS.getSegment error: incorrect input args\n');
				SegmentInfo = [];
				SourceInfo = [];
				return
			end
			% if segList is empty, use the full SegmentList
			if isempty(segList)
				segList = obj.SegmentList;
			end
			Nseg = length(segList);
			%----------------------------------------------------
			% open file
			%----------------------------------------------------
			obj.openFile;
			%----------------------------------------------------
			% pre-allocate things
			%----------------------------------------------------
			t = zeros(Nseg, 1);
			%----------------------------------------------------
			% loop through segments
			%----------------------------------------------------
			for n = 1:Nseg
				% get info and source info
				[t(n), SegmentInfo(n)] = ...
						ns_GetSegmentInfo(obj.Hfile, segList(n)); %#ok<AGROW>
				[t(n), SourceInfo(n)] = ...
						ns_GetSegmentSourceInfo(obj.Hfile, segList(n), 1); %#ok<AGROW>
			end
			%----------------------------------------------------
			% close file
			%----------------------------------------------------
			obj.closeFile;
		end
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		
		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function N = getNeural(obj)
		%------------------------------------------------------------------------
		% Neural = NS.getNeural
		%------------------------------------------------------------------------
		% read neural entities
		%------------------------------------------------------------------------
		% Neural Event Entities (4)
		%
		% Timestamps of event and segment entitities that
		% are known to represent neural action potential firing times.  For
		% example, if a segment entity contains sorted neural spike waveforms, 
		% each sorted unit is also exported as a neural entity.
		% Each index of a neural event entity refers to a timestamp for each
		% neural event.  Each neural event entity contains event times for a
		% single neural source.  The number of indices is equal to the number
		% of entries for that neural event entity.
		%------------------------------------------------------------------------

			% checks
			N = [];
			if obj.checkReadStatus == 0
				return
			end
			if (obj.nNeural == 0)
				disp('No Neural entities available!');
				return
			end
			% preallocate N
			N = repmat(	struct(	'Info', [], ...
										'Data', [], ...
										'TimeStamps', [], ...
										'ItemCount', [] ), ...
							obj.nNeural, 1);
			s = zeros(obj.nNeural, 1);
			% built list of labels
			NeuralLabels = char(obj.EntityInfo(obj.NeuralList).EntityLabel);
			% loop through neural list
			for n = 1:obj.nNeural
				% Retrieve the data
				[s(n), N(n).Info] = ns_GetNeuralInfo(obj.Hfile, obj.NeuralList(n));
				[s(n), N(n).Data] = ...
					ns_GetNeuralData(obj.Hfile, ...
											obj.NeuralList(n), ...
											1, ...
											obj.EntityInfo(obj.NeuralList(n)).ItemCount);

				indx = 1:(size(N(n).Data, 1) * size(N(n).Data, 2));

				% Get the neural timestamps
				N(n).Data = reshape(N(n).Data, size(N(n).Data, 1) * size(N(n).Data, 2), 1);
				N(n).TimeStamps = {N(n).Data(indx)};

				% Match the neural events with their unit ID
				% Have to figure out which Neural entities correspond 
				% with the selected segment entities
				%list = strmatch(EntityInfo(SegmentList(channel(cChannel))).EntityLabel, NeuralLabels, 'exact');
				%tempID = ones(length(indx) / length(list), 1) * [tmp.Info(:).SourceUnitID];
				%N(n).Units = {tempID(:)};

				% Remember how many neural events were found
				N(n).ItemCount = length(N(n).TimeStamps);
			end
		end	% END getNeural
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		
	end	% END of PUBLIC METHODS
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	
	
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% PROTECTED methods
	%------------------------------------------------------------------------
	methods (Access = protected)
	
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function obj = openFile(obj)
		%------------------------------------------------------------------------
		% NS.openFile(obj)
		%------------------------------------------------------------------------
		% Open data file and get some info about the file
		%------------------------------------------------------------------------

			%--------------------------------
			% perform some checks
			%--------------------------------
			if ~isempty(obj.Hfile)
				% file is already open!
				return
			end
			if isempty(obj.datafullname)
				% no data file set
				fprintf('NS error: data file name not set!\n')
				return
			end
			if ~obj.libraryLoaded
				% Library not loaded - load it now
				fprintf('NS warning: library not loaded!\n')
				fprintf('\t loading default library %s\n', obj.DLLfullname);
				obj.setLibrary;
			end			
			%---------------------------------------
			% open file using NeuroShare function
			%---------------------------------------
			[nsresult, obj.Hfile] = ns_OpenFile(obj.datafullname);
			% check 
			if (nsresult ~= 0)
				fprintf('%s: Data file %s did not open!\n', ...
							mfilename, obj.datafullname);
				obj.Hfile = 0;
				return
			end
		end	% END openFile
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function obj = closeFile(obj)
		%------------------------------------------------------------------------
		% NS.closeFile
		%------------------------------------------------------------------------
		% Close data file
		%------------------------------------------------------------------------
			if isempty(obj.Hfile)
				fprintf('File interface already closed!\n')
				return
			else
				ns_CloseFile(obj.Hfile);
				obj.Hfile = [];
			end
		end	% END closeFile
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function status = checkReadStatus(obj)
		%------------------------------------------------------------------------
		% status = NS.checkReadStatus
		%------------------------------------------------------------------------
		% use to check IO status
		%------------------------------------------------------------------------
			if ~obj.libraryLoaded
				fprintf('NS: library not loaded!\n');
				status = 0;
				return
			elseif isempty(obj.Hfile)
				fprintf('NS: file not open\n');
				status = 0;
				return
			else
				status = 1;
			end
		end
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------		
		
	end	% END of PROTECTED METHODS
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	
	
	
	
	
	
	
	
	
	
end	% END of NS
%*****************************************************************************
%*****************************************************************************
%*****************************************************************************
