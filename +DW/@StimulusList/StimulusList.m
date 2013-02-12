%-----------------------------------------------------------------------------
% StimulusList.m
%-----------------------------------------------------------------------------
% DataMat Toolbox
% DW package
% Class Definition
%-----------------------------------------------------------------------------
%	Properties:
% 		N						# of Different stimuli
% 		MarkerList			{N X 1} list of marker indices
% 		Tagstring			{N X 1} list of unique marker strings
% 		S = {};				{N X # audio chan} list of Stimulus objects (by type)
% 		Type					{N X 2} list of Stimulus types ('WAV', 'TONE', etc)
% 		Channel				[N X 1] character array of channel ('0' 'L' 'R' 'B')
% 		Var
% 		Nsweeps				[N X 1] array of sweeps/stim
% 		Sweepstart			{N X 1} list of vectors of sweep start times in usec
%									each vector will have a time for each type
% 									and sweep
% 		Sweepend				same format as Sweepstart, but contains sweep end 
% 									times in usec
% 		PreSweepMark		{N X 1} list of indices into Markers (see DW.Data)
% 									for sweep before 
% 		PostSweepMark		{N X 1} list of indices into Markers for sweep after
%-----------------------------------------------------------------------------
% See also: DWdata, Marker, Probe, Unit
%-----------------------------------------------------------------------------

%-----------------------------------------------------------------------------
%	Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%-----------------------------------------------------------------------------
% Created: 18 January, 2012 (SJS)
%
% Revisions:
%-----------------------------------------------------------------------------
% TO DO:
%-----------------------------------------------------------------------------

%*****************************************************************************
%*****************************************************************************
%*****************************************************************************
% class definition
%*****************************************************************************
%*****************************************************************************
%*****************************************************************************
% classdef (ConstructOnLoad = true) StimulusList < handle
classdef StimulusList < handle
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Properties
	%------------------------------------------------------------------------
	properties
		N
		MarkerList
		Tagstring
		S = {};
		GroupList = {};
		Type
		Channel
		Var
		Nsweeps
		Sweepstart
		Sweepend
		PreSweepMark
		PostSweepMark
	end	% end of properties
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------

	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Define methods
	%------------------------------------------------------------------------
	methods		
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		% Constructor Method
		%---------------------------------------------------------------------
		function obj = StimulusList(varargin)
		%---------------------------------------------------------------------	
		% StimulusList
		% Object Constructor method
		%---------------------------------------------------------------------
		% StimulusList()	when called with no arguments, returns empty
		%						StimulusList object.  
		% SimulusList(Markers)	When given a list of marker objects, the 
		%								markers will be analyzed for specific
		%								stimulus types and then sorted.
		%---------------------------------------------------------------------
			if nargin == 0
				% if no arguments were provided, return empty StimulusList
				fprintf('%s: building empty StimulusList\n', mfilename);
				return
			else
				fprintf('%s: building StimulusList\n', mfilename);
				% if Marker objects provided, parse them to get stimulus info
				obj.parseMarkersIntoStimuli(varargin{1});
			end
		end	% END StimulusList CONSTRUCTOR
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		% General Methods
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		function parseMarkersIntoStimuli(obj, Markers)
		%---------------------------------------------------------------------
		% StimulusList.parseMarkersIntoStimuli(Markers)
		%---------------------------------------------------------------------
			DW.DataWaveDefaults;	% load defaults
			%-----------------------------------------------------------
			% find unique markers (stimuli)
			%-----------------------------------------------------------
			fprintf('%s: finding unique markers...\n', mfilename);
			% # of markers
			Nmarkers = length(Markers);
			% need to get list of all stimuli in text form
			% preallocate stim tag list
			stimtext = cell(Nmarkers, 1);
			% loop through markers
			for n = 1:Nmarkers
				% convert each marker's tag values to a string...
				tmp = cell2str(Markers(n).getStimulus');
				% ...and store in stimtext 
				stimtext{n} = tmp{1};
			end	% END n
			% then find unique values of strings, store 
			% Tagstring, Indices, Nstimuli
			[obj.Tagstring, obj.MarkerList, obj.N] = findUniqueText(stimtext);
			obj.Tagstring = obj.Tagstring';
			clear stimtext;
			fprintf('\t found %d unique markers from %d total markers\n', ...
									obj.N, Nmarkers);

			%-----------------------------------------------------------
			% determine type and channel from markers
			%-----------------------------------------------------------
			fprintf('%s: finding stimulus types and channels...\n', mfilename);
			% create, and assign type
			obj.Type = cell(obj.N, 2);
			obj.Channel = char(obj.N, 1);
			for n = 1:obj.N
				% check stimuli by examining the first element in each
				% marker/stimulus subtype (this is assuming that the sorting of
				% the markers went as planned!!!)
				rstim = Markers(obj.MarkerList{n}(1)).SoundTypeR;
				lstim = Markers(obj.MarkerList{n}(1)).SoundTypeL;
				% check if either or both of rstim == 0
				if strcmpi(rstim, '0') && strcmpi(lstim, '0')
					% if both == 0, set Channel to '0' (none)
					obj.Channel(n) = '0';
				elseif strcmpi(rstim, '0')
					% if rstim == 0, then Channel must be 'L'
					obj.Channel(n) = 'L';
				elseif strcmpi(lstim, '0')
					% if lstim == 0, then Channel must be 'R'
					obj.Channel(n) = 'R';
				else
					% otherwise, both channels active, Channel set to 'B'
					obj.Channel(n) = 'B';
				end
				fprintf('stim %d: channel = %s', n, obj.Channel(n));				
				fprintf('\t\t(%s\t%s)\n', lstim, rstim);
				% set Type
				obj.Type(n, :) = {lstim rstim};
			end	% END n
			
			%-----------------------------------------------------------
			% sort into Stimuli
			%-----------------------------------------------------------
			fprintf('%s: sorting into stimuli...\n', mfilename);
			% allocate S(timulus) cell array to hold object list
			obj.S = cell(obj.N, 2);
			% determine stimulus type and then initialize/assign object
			for n = 1:obj.N
				for c = L:R
					% assign channel object
					switch upper(obj.Type{n, c})
						case '0'		% no sound
							obj.S{n, c} = [];
						case 'BBN'	% noise
							obj.S{n, c} = ...
												DW.Noise(Markers(obj.MarkerList{n}(1)));
						case {'WAVFILE', 'WAV'}	% wav file
							obj.S{n, c} = ...
												DW.Wav(Markers(obj.MarkerList{n}(1)));
						case 'TONE'	% tone
							obj.S{n, c} = ...
												DW.Tone(Markers(obj.MarkerList{n}(1)));
						otherwise
							warning('parseMarkersIntoStimuli: unknown stim type %s', ...
												mfilename, obj.Type{n}{c});
					end	% END switch
				end	% END c
			end	% END n
			fprintf('\nMarker Parsing complete\n');
			
			% count # of sweeps/stim
			obj.Nsweeps = zeros(obj.N, 1);
			for n = 1:obj.N
				obj.Nsweeps(n) = length(obj.MarkerList{n});
			end
		end		% END parseMarkersIntoStimuli 
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------

		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		function extractTimeFromMarkers(obj, Markers, varargin)
		%---------------------------------------------------------------------
		% StimulusList.extractTimeFromMarkers(Markers)
		% StimulusList.extractTimeFromMarkers(Markers, 'SweepDuration', SweepDur)
		% StimulusList.extractTimeFromMarkers(Markers, 'MaxTime', MaxTime)
		%---------------------------------------------------------------------
			DW.DataWaveDefaults;	% load defaults
			%-----------------------------------------------------------
			%-----------------------------------------------------------
			SweepDuration = [];
			MaxTime = [];
			if ~isempty(varargin)
				vindx = 1;
				while vindx <= length(varargin)
					switch upper(varargin{vindx})
						case 'SWEEPDURATION'
							SweepDuration = varargin{vindx+1};
							vindx = vindx + 2;
						case 'MAXTIME'
							MaxTime =  varargin{vindx+1};
							vindx = vindx + 2;
						otherwise
							error('%s: unknown option %s', varargin{vindx});
					end
				end
			end
			if isempty(SweepDuration)
				SweepDuration = SWEEP_DURATION;
			end
			
			%-----------------------------------------------------------
			% init
			%-----------------------------------------------------------
			fprintf('%s: calculating sweep times...\n', mfilename);
			% # of markers
			Nmarkers = length(Markers);
			% allocate storage
			obj.Sweepstart = cell(obj.N, 1);
			obj.Sweepend = cell(obj.N, 1);
			obj.PreSweepMark = cell(obj.N, 1);
			obj.PostSweepMark = cell(obj.N, 1);

			for s = 1:obj.N
				fprintf('%d sweeps for stimulus %d\n', obj.Nsweeps(s), s);
				obj.Sweepstart{s} = zeros(obj.Nsweeps(s), 1);
				obj.Sweepend{s} = zeros(obj.Nsweeps(s), 1);
				obj.PreSweepMark{s} = zeros(obj.Nsweeps(s), 1);
				obj.PostSweepMark{s} = zeros(obj.Nsweeps(s), 1);

				for w = 1:obj.Nsweeps(s)
					m = obj.MarkerList{s}(w);
					fprintf('\t%.0f\t\t%.0f\t%s\t%.0f\n', ...
						Markers(m).OutputTimestampL, ...
						Markers(m).OutputTimestampR, Markers(m).WavFilenameR, ...
						Markers(m).AttenuationR);

					if m == 1
						preindx = 0;
					else
						preindx = m - 1;
					end
					if m == Nmarkers
						postindx = Nmarkers;
					else
						postindx = m + 1;
					end

					if obj.Channel(s) == 'L'
						obj.Sweepstart{s}(w) = Markers(m).OutputTimestampL;
						if isempty(SweepDuration)
							obj.Sweepend{s}(w) = ...
															Markers(postindx).OutputTimestampL;
						else
							obj.Sweepend{s}(w) = obj.Sweepstart{s}(w) + SweepDuration;
						end
					elseif obj.Channel(s) == 'R'
						obj.Sweepstart{s}(w) = Markers(m).OutputTimestampR;
						if isempty(SweepDuration)
							obj.Sweepend{s}(w) = Markers(postindx).OutputTimestampR;
						else
							obj.Sweepend{s}(w) = obj.Sweepstart{s}(w) + SweepDuration;
						end
					elseif obj.Channel(s) == 'B'
						% use min value of l and r
						obj.Sweepstart{s}(w) = ...
													min([Markers(m).OutputTimestampL ...
															Markers(m).OutputTimestampR]);
						if isempty(SweepDuration)
							obj.Sweepend{s}(w) = ...
												min([Markers(postindx).OutputTimestampL ...
														Markers(postindx).OutputTimestampR]);
						else
							obj.Sweepend{s}(w) = obj.Sweepstart{s}(w) + SweepDuration;
						end
					else
						error('%s: bad channel value %s', mfilename, obj.Channel(s));
					end
					% correct final sweep time
					if postindx == Nmarkers
						if isempty(MaxTime)
							obj.Sweepend{s}(w) = ...
									obj.Sweepstart{s}(w) + SweepDuration;
						else
							obj.Sweepend{s}(w) = MaxTime;
						end
					end
					obj.PreSweepMark{s}(w) = preindx;
					obj.PostSweepMark{s}(w) = postindx;
				end	% END sweep
			end	% END stim

		end	% END extractTimeFromMarkers FUNCTION
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		
		%---------------------------------------------------------------------
		function varargout = findCommon(obj)
		%---------------------------------------------------------------------
		% varargout = StimulusList.findCommon
		%---------------------------------------------------------------------
			
			%----------------------------------
			% initialize values
			%----------------------------------
			obj.GroupList = {};
			runFlag = 1;
			nUnique = 0;
			% masterList is list of indices into rows of obj.S
			masterList = 1:obj.N;
			% search all indices initially - this will be shrunk as
			% matches are found
			searchList = 1:obj.N;
			%----------------------------------
			% loop
			%----------------------------------
			while runFlag
				% use first element in searchList as test
				stim = searchList(1);
				% check channel for current stimulus
				c = obj.Channel(stim);
				% compare current stim to other stims
				if strcmpi(c, 'L') || strcmpi(c, 'B')
					% if channel is L or Both, compare left, store matches in lcomp
					[lcomp, llist] = obj.S{stim, 1}.match(obj.S(searchList, 1));
				end
				if strcmpi(c, 'R') || strcmpi(c, 'B')
					% if channel is R or Both, compare right, store matches in rcomp
					[rcomp, rlist] = obj.S{stim, 2}.match(obj.S(searchList, 2));
				end
				if strcmpi(c, 'B')
					% if channel is Both, AND the lcomp and rcomp
					comp = lcomp & rcomp;
				else
					% otherwise, use appropriate channels results
					if strcmpi(c, 'L')
						comp = logical(lcomp);
					else
						comp = logical(rcomp);
					end
				end
				% store unique indices (from masterList) in GroupList
				nUnique = nUnique + 1;
				obj.GroupList{nUnique} = masterList(comp);
				% eliminate them from the list to search
				searchList = searchList(~(comp));
				% and shrink the unsearched masterList as well
				masterList = masterList(~(comp));
				% Check if we're done
				if isempty(searchList)
					runFlag = 0;
				end
			end
			if nargout
				varargout{1} = obj.GroupList;
			end
		end	% END findMatchingStim FUNCTION
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------

		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		% Overloaded Methods
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------

		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		% set/get Methods
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		
	
	end	% End of methods
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
end	% End of classdef
%*****************************************************************************
%*****************************************************************************
%*****************************************************************************




