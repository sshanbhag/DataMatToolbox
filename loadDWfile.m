function varargout = loadDWfile(varargin)
%------------------------------------------------------------------------
% [D, Stimulus, errFlg, rawdata] = loadDWfile(fname, <options>)
%------------------------------------------------------------------------
% DataMat toolbox
%------------------------------------------------------------------------
% loads and processes data stored in Wenstrup lab Datawave .txt files
% 
%------------------------------------------------------------------------
% Input Arguments:
% 	fname		path to and name of DataWave exported text file 
% 				(usually with .dat file extension)
%				If fname is not provided, a GUI window will open to select file.
% 
% 	Options:		given as a pair of '<option name>', '<option value>'
% 	
% 		'SaveMatfile'		'Yes' -> Saves data in .mat file with same name as
% 											data file
% 								'No'	-> does not save data, does not query user
%								'<filename>' -> mat file will be written to filename.mat
% 								
% 								If option is not provided as input, the 
% 								program will ask the user is data should be saved 
% 								to a .mat file
%
%		'OutputPath'		<path> -> path to store output mat file
% 								
% 		'PlotData'			'Yes' -> Plots of data will be generated
% 								'No'	-> Plots of data will not be displayed, 
% 								         and user will not be queried
% 											
% 								If option is not provided, the program will ask if
% 								plots should be displayed
% 								
% Output Arguments:
% 	D			output structure
%
% 		info				DataWave file information structure:
%			file				filename
%			path				path to file = pname;
%			Nlines			# of lines (including header) in file
%			header			header structure
%									line			cell vector of header line text
%									fields		cell vector of field names
%									nfields		# of tab-delimited fields in each header line
%			data1				text from data line 1
%			ndata1			# of fields in data line 1
% 
% 		Marker			DataWave event marker structure, containing the following
%							arrays [Nmarkers X 1]:
% 			string						raw string from text file
% 			Timestamp					time of event, in microseconds
% 			id								?
% 			OutputFilename				output data file name
% 			AttenuationR				R channel attenuation setting (dB)
% 			BBNlowerFreqR				if noise stimulus, 
% 			BBNupperFreqR
% 			AmplitudeR
% 			TimeShiftR
% 			RampUpR
% 			HoldTimeR
% 			RampDownR
% 			OutputTimestampR
% 			OutputTimeWithDelayR
% 			FixedDelayR
% 			PA5idR
% 			WavFilenameR
% 			ToneFreqR
% 			PhaseDegR
% 			AttenuationL
% 			BBNlowerFreqL
% 			BBNupperFreqL
% 			AmplitudeL
% 			TimeShiftL
% 			RampUpL
% 			HoldTimeL
% 			RampDownL
% 			OutputTimestampL
% 			OutputTimeWithDelayL
% 			FixedDelayL
% 			PA5idL
% 			WavFilenameL
% 			ToneFreqL
% 			PhaseDegL
% 			M
% 			StimulusTypeR
% 			StimulusTypeL
% 			wavFilesR
% 			wavFilesL
% 		
%		MarkerTimes		[# of markers X 1] vector of Marker timestamps, in
%							microseconds
% 		
% 		Probe				spike information structure, cluster information added
% 			t					timestamp vector (microseconds)
% 			cluster			cluster ID number for spikes 
% 									0 => unsorted spike
% 									1, 2, ... N => sorted spike cluster number
% 			nclusters		number of clusters on this probe (electrode or tetrode)
% 
%		UnitData			UnitData struct array holds spiketimes, probe information 
% 							and unit ID for sorted (or unsorted) spikes
% 			probe				probe ID #
% 			unit				unit ID # (usually equal to n-1)
% 			indices			vector of indexes into master event list [1xN double]
% 			timestamp		vector of time stamps, microseconds [1xN double]
% 			sorted			1 if spikes are sorted/thresholded, 0 for unsorted spikes
% 
% 	Stimulus			Structure array, with following elements
% 		Type							'TONE', 'NOISE', 'WAVFILE'
% 		Channel						'R', 'L', 'B'
% 		Indices						array of indices into Marker arrays (above)
% 		Nreps							# of times this stimulus was presented 
% 		Var							varying variables for this stimulus
% 		Timestamp					first occurance of this Stimulus
% 		id
% 		OutputFilename
% 		AttenuationR
% 		BBNlowerFreqR
% 		BBNupperFreqR
% 		AmplitudeR
% 		TimeShiftR
% 		RampUpR
% 		HoldTimeR
% 		RampDownR
% 		OutputTimestampR
% 		OutputTimeWithDelayR
% 		FixedDelayR
% 		PA5idR
% 		WavFilenameR
% 		ToneFreqR
% 		PhaseDegR
% 		AttenuationL
% 		BBNlowerFreqL
% 		BBNupperFreqL
% 		AmplitudeL
% 		TimeShiftL
% 		RampUpL
% 		HoldTimeL
% 		RampDownL
% 		OutputTimestampL
% 		OutputTimeWithDelayL
% 		FixedDelayL
% 		PA5idL
% 		WavFilenameL
% 		ToneFreqL
% 		PhaseDegL
% 		Tagstring
% 		Nsweeps						# of times this stimulus was presented
% 		Sweepstart					list of start times for this stimulus
% 		Sweepend						list of end times for this stimulus
% 		LAttenVals					Left channel attenuation values
% 		LAttenIndices				Indices for respective attenuation sweeps
% 		RAttenVals					Right channel attenuation values
% 		RAttenIndices				indices for R attenuaton sweeps
% 		Spiketimes					{# units, Nsweeps} cell array of unit spike
%										timestamps
% 
%	errFlg			Error flag
% 			0					no error
% 			1					user cancelled file opening
% 			2					no fields found in header lines
% 			3					file not found
%
%	rawdata			"raw" data cell array
%------------------------------------------------------------------------
% See also: readDWfileinfo 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 24 January, 2011 (SJS)
%
% Revisions:
%	17 Feb 2011 (SJS)
% 	 - finished parsing of data{} into Probe and Marker structs
% 	 -	converted to function
% 	 -	documentation
%	27 Apr 2011 (SJS): added some user-feedback bits
%	13 May 2011 (SJS): fixed problem with use of fully-spec'ed filenames
%  25 May, 2011 (SJS): adapting for new data format
%	5 July, 2011 (SJS): documenting and commenting
%	28 July, 2011 (SJS): added Hout to outputs
%	4 August, 2011 (SJS): 
% 		-	added Background retrieval and item to outputs
% 		-	removed path name as input
%	26 January, 2012 (SJS): added 'OutputPath' option
%	17 May, 2012 (SJS): working on objectifying
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

%-----------------------------------------------------------
% load defaults
%-----------------------------------------------------------
DataWaveDefaults;

%-----------------------------------------------------------
% check input arguments, act depending on inputs
%-----------------------------------------------------------
if nargin == 0
	% no filename or path provided
	% open ui panel to get filename and path from user
	[fname, pname] = uigetfile('*.txt', 'Select Datawave Exported Text File');
	% return if user pressed cancel button
	if isequal(fname, 0) || isequal(pname,0)
		disp('Cancelled...')
		for n = 1:nargout
			varargout{n} = [];
		end
		if nargout >= 4
			varargout{4} = 1;
		end
		return
	end
end

if nargin > 0
	% only filename was provided, assume path is included or that file is
	% in current directory 
	[pname, fname, ext] = fileparts(varargin{1});
	if isempty(pname)
		pname = pwd;
	end
	fname = [fname ext];
end


if nargin >= 2
	chkargs = varargin(2:end);
	for n = 1:length(chkargs)
		if ischar(chkargs{n})
			if strcmpi(chkargs{n}, 'SaveMatfile')
				if strncmpi(chkargs{n+1}, 'Y', 1)
					saveMatFlg = 1;
				elseif strncmpi(chkargs{n+1}, 'N', 1)
					saveMatFlg = 0;
				else
					saveMatFlg = 1;
					matfname = chkargs{n+1};
				end
			elseif strcmpi(chkargs{n}, 'PlotData')
				if strncmpi(chkargs{n+1}, 'Y', 1)
					plotDataFlg = 1;
				else
					plotDataFlg = 0;
				end
			elseif strcmpi(chkargs{n}, 'OutputPath')
				if ischar(chkargs{n+1})
					OutputPath = chkargs{n+1};
				end				
			end
		end
	end
end

%-----------------------------------------------------------
% get file info, read header
%-----------------------------------------------------------
disp(['Reading information from file ' fname ' ...']);
[dwinfo, errFlg] = readDataWaveTextInfo(fname, pname);

% check if errFlg ~= 0 (error detected)
if errFlg
	warning('DWFILE:ReadError', '%s readDataWaveTextInfo returned error flag %d', ...
												mfilename, errFlg);
end

%-----------------------------------------------------------
% perform some checks on file information, if okay, read in
% data from DW text file
%-----------------------------------------------------------
if ~dwinfo.header.nfields(1)
	save('loadDWfile.error.mat', 'dwinfo', '-MAT');
	error('%s: no fields found in header line 2 of data file', mfilename);
elseif ~dwinfo.Ncols
	save('loadDWfile.error.mat', 'dwinfo', '-MAT');
	error('%s: no fields found in data line 1 of data file', mfilename);
end

%-----------------------------------------------------------
% parse header data to get information about data set
%-----------------------------------------------------------
disp('Parsing header information...');
[dwinfo, errFlg] = parseDataWaveTextHeader(dwinfo);

%-----------------------------------------------------------
% read in data from file
%-----------------------------------------------------------
disp(['Reading data from file ' dwinfo.filename ' ...']);
[rawdata, errFlg] = readDataWaveTextFile(dwinfo);

%-----------------------------------------------------------
% get Marker information and retrieve markers from data
%-----------------------------------------------------------
disp( 'Processing Marker information ...')

% check for markers
if ~dwinfo.NMarkerCols
	% if no marker cols, error!
	error('%s: no markers detected in file', mfilename)
end

% Pull in Marker data
disp('Reading Marker Data...')

% first and last columns of Marker information
mStartCol = dwinfo.MarkerCols(1);
mEndCol = dwinfo.NMarkerCols;

% initialize marker Counter variable
markerCount = 0;

% loop, using # of data lines as upper bound
for L = 1:dwinfo.Ndatalines
	% check if marker information is present
	if isempty(rawdata{L}{mStartCol})
		% if not, break
		disp(['found end of Marker Data, line ' num2str(L)])
		break;
	else
		% store marker information
		
		% increment marker index
		markerCount = markerCount + 1;

		% save time stamp as number
		MarkerData(markerCount).t = str2double(rawdata{L}{mStartCol});
		% save times in stand-alone vector
		MarkerTimes(markerCount) = MarkerData(markerCount).t;
		
		% initialize fields for marker data
		% I am assuming two potential types of data
		%	text		some sort of string, e.g. .wav file name
		%	value		a numeric value (e.g., attenuation value)
		MarkerData(markerCount).string = rawdata{L}(mStartCol:mEndCol);
	end
end

% Toss out late Marker Line as it is only a place holder
%	(added 1/26/2012, SJS)
MarkerData = MarkerData(1:(markerCount - 1));
MarkerTimes = MarkerTimes(1:(markerCount - 1));
markerCount = markerCount - 1;

% process marker information into arrays
[Marker, errFlg] = parseDataWaveMarkers(MarkerData, dwinfo);
if errFlg
	error('%s: parseDataWaveMarkers error # %d', mfilename, errFlg)
end

%-----------------------------------------------------------
% get # of markers by looking at length of the Timestamp vector 
% could use any of the vectors
%-----------------------------------------------------------
Marker.Nmarkers = length(Marker.Timestamp);

%-----------------------------------------------------------
% Pull in Spike Channel Data
%-----------------------------------------------------------
% check # of Spike channels
if ~dwinfo.NSpikeCols
	% if 0, error
	error('%s: no spike data channels detected in header', mfilename)
else
	% otherwise, build Probe data structure
	for n = 1:dwinfo.NSpikeCols
		ProbeData(n) = struct('t', [], 'cluster', []);
	end
end

disp('Reading and Parsing Probe Data...')
% loop through spike columns (i.e., tetrodes)
for p = 1:dwinfo.NSpikeCols
	% get current column number for spike data
	c = dwinfo.SpikeCols(p);
	% pull out the data for this column
	l = 0;
	tmpline = rawdata{1};
	loopFlag = 1;
	while ~isempty(tmpline{c}) && (l < dwinfo.Ndatalines) && loopFlag
		l = l+1;
		tmpline = rawdata{l};
		if ~isempty(tmpline{c})
			ProbeData(p).t(l) = str2double(tmpline{c});
			ProbeData(p).cluster(l) = str2double(tmpline{c+1});
		else
			loopFlag = 0;
		end
	end
end

%-----------------------------------------------------------
% Parse ProbeData, organizing by Probe and unit
%-----------------------------------------------------------
disp('Parsing ProbeData...')
[UnitData, ProbeData, errFlag] = parseDataWaveProbes(ProbeData, Marker);

dwinfo.Nunits = length(UnitData);

%-----------------------------------------------------------
% Create output structure
%-----------------------------------------------------------
D = struct(	'Info', dwinfo, ...
				'Probe', ProbeData, ...
				'Marker', Marker, ...
				'MarkerTimes', MarkerTimes, ...
				'UnitData', UnitData ...
			);
		
%-----------------------------------------------------------
% Create Stimulus structure
%-----------------------------------------------------------
disp('Creating Stimulus structure array...');
Stimulus = buildStimulusStruct(D);

%-----------------------------------------------------------
% Create Background structure
%-----------------------------------------------------------
disp('Creating Background structure array...');
Background = buildBackgroundStruct(D);

%-----------------------------------------------------------
% Create UnitInfo struct and add to Info
%-----------------------------------------------------------
[D.Info.UnitInfo, D.Info.UnitList] = buildUnitInfo(D, Stimulus);

%-----------------------------------------------------------
% save data to mat file
%-----------------------------------------------------------
% see if matfname is defined (i.e., was provided in function call)
if ~exist('matfname', 'var')
	% if not, create it from the .txt file name
	[~, tmpfname, ext] = fileparts(fname);
	% remove the _Sheetmaker from the filename
	tmpfname = regexprep(tmpfname, '(_Sheetmaker)', '');
	% remove the _spksorted from the filename
	tmpfname = regexprep(tmpfname, '(_spksorted)', '');
	% append .mat
	matfname = [tmpfname '.mat'];
	clear tmpfname
end

% check if path was provided
if exist('OutputPath', 'var')
	if ~exist(OutputPath, 'dir')
		% create output path if it doesn't exist
		warning('%s: creating output directory %s', mfilename, OutputPath);
		mkdir(OutputPath);
	end
else
	OutputPath = [];
end

% prepend matfname with OutputPath
matfname = fullfile(OutputPath, matfname);

% if saveMatFlg is set, force write of matfile
if exist('saveMatFlg', 'var')
	if saveMatFlg
		save(matfname, 'D', 'Stimulus', 'Background', '-MAT');
	end
else
	% user didn't specify preference on data file saving in input args
	% check if file exists
	if exist(matfname, 'file')
		disp([mfilename ':  warning, output .mat file ' matfname ' alread exists!!!']);
		q = input('Overwrite file? (y/n) ', 's');
		if ~isempty(q)
			if strcmpi(q(1), 'Y')
				save(matfname, 'D', 'Stimulus', 'Background', '-MAT')
			end
		end
	else
		save(matfname, 'D', 'Stimulus', 'Background', '-MAT')
	end
end

%-----------------------------------------------------------
% assign outputs depending on number of outputs requested
%-----------------------------------------------------------
if any(nargout == [0 1 2 3 4 5])
	varargout{1} = D;
end

if any(nargout == [2 3 4 5])
	varargout{2} = Stimulus;
end

if any(nargout == [3 4 5])
	varargout{3} = Background;
end

if any(nargout == [4 5])
	varargout{4} = errFlg;
end

if nargout == 5
	varargout{5} = rawdata;
end


if exist('plotDataFlg', 'var')
	if plotDataFlg
		if ~isempty(regexp(D.Info.file, '(FRA)', 'once'))
			plotFRAData
		else
			Hout = plotData(D, Stimulus);
		end
	else
		Hout = [];
	end
else
	% user didn't specify preference on data plotting in input args
	q = input('Plot rasters and paths (y/n)? ', 's');
	if ~isempty(q)
		if strcmpi(q(1), 'Y')
			disp('Plotting rasters and psths for all units...')

			if ~isempty(regexp(D.Info.file, '(FRA)', 'once'))
				plotFRAData
			else
				Hout = plotData(D, Stimulus);
			end
		end
	else
		Hout = [];
	end
end

% return plot handles if requested
if nargout == 6
	varargout{6} = Hout;
end


