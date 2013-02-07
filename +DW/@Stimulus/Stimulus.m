%-----------------------------------------------------------------------------
% Stimulus.m
%-----------------------------------------------------------------------------
% DataMat Toolbox
% DW package
% Class Definition
%-----------------------------------------------------------------------------
%	Properties:
% 		Type							'TONE', 'NOISE', 'WAVFILE'
% 		Channel						'R', 'L', 'B'
% 		Indices						array of indices into Marker arrays (above)
% 		Nreps							# of times this stimulus was presented 
% 		Var							varying variables for this stimulus
% 		Nsweeps						# of times this stimulus was presented
% 		Sweepstart					list of start times for this stimulus
% 		Sweepend						list of end times for this stimulus
% 		LAttenVals					Left channel attenuation values
% 		LAttenIndices				Indices for respective attenuation sweeps
% 		RAttenVals					Right channel attenuation values
% 		RAttenIndices				indices for R attenuaton sweeps
% 		Spiketimes					{# units, Nsweeps} cell array of unit spike
% 										timestamps
% 		Timestamp					first occurance of this Stimulus
%		Tagstring					for debugging, will disappear eventually
%-----------------------------------------------------------------------------
% See also: DWdata, Marker, Probe, Unit
%-----------------------------------------------------------------------------

%-----------------------------------------------------------------------------
%	Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%-----------------------------------------------------------------------------
% Created: 4 June, 2012 (SJS)
%
% Revisions:
%-----------------------------------------------------------------------------
% TO DO:
%
%-----------------------------------------------------------------------------

%*****************************************************************************
%*****************************************************************************
%*****************************************************************************
% class definition
%*****************************************************************************
%*****************************************************************************
%*****************************************************************************
% classdef (ConstructOnLoad = true) Stimulus < handle
classdef Stimulus < handle
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Protected Properties
	%------------------------------------------------------------------------
	properties (SetAccess = protected)
		Type
		Channel
		Amplitude
		Attenuation
		TimeShift
		RampUp
		HoldTime
		RampDown
		OutputTimestamp
		OutputTimeWithDelay
		FixedDelay
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
		function obj = Stimulus(varargin)
		%---------------------------------------------------------------------	
		% Stimulus
		% Constructor method
		%---------------------------------------------------------------------
		% Stimulus()	when called with no arguments, returns empty
		%					Stimulus object
		%---------------------------------------------------------------------
			DW.DataWaveDefaults;
			%--------------------------------------------------------
			% Check inputs
			%--------------------------------------------------------
			% if no inputs, return
			if isempty(varargin)
				return
			end
			% if 1 input, need to perform a few checks
			if length(varargin) == 1
				% check if varargin{1} was directly passed in from
				% a subclass (using obj@DW.Stimulus(varargin) syntax)
				if strcmpi(inputname(1), 'varargin')
					% make local copy of varargin{1}
					tmp = varargin{1};
					if length(tmp) == 1
						% tmp is marker only, assume RIGHT channel
						obj.buildStimulusFromMarker(tmp{1}, R);
					elseif length(tmp) == 2
						% caller provided Marker and Channel as inputs
						obj.buildStimulusFromMarker(tmp{1}, tmp{2});
					else
						fprintf('%s: strange inputs... ', mfilename);
						fprintf('%s\t', tmp);
						fprintf('\n');
						error('%s: invalid inputs', mfilename);
					end
					clear tmp
				else
					% if no channel provided, assume channel it RIGHT
					obj.buildStimulusFromMarker(varargin{1}, R);
				end
			elseif length(varargin) == 2
				% caller provided Marker and Channel as inputs
				obj.buildStimulusFromMarker(varargin{1}, varargin{2});
			else
				error('%s: invalid input args', mfilename);
			end

		end		% END Stimulus constructor
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		% General Methods
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		function obj = buildStimulusFromMarker(obj, M, C)
		%---------------------------------------------------------------------
		% Stimulus.buildStimulusFromMarker(M, C)
		%---------------------------------------------------------------------
		% given marker object M and channel C (L == 1, R == 2), assign values
		% from marker to appropriate values in Stimulus object
		%---------------------------------------------------------------------
			DW.DataWaveDefaults;
			markerbase = {	'Amplitude', 'Attenuation', 'TimeShift', ...
							'RampUp', 'HoldTime', 'RampDown', ...
							'FixedDelay', 'OutputTimestamp', 'OutputTimeWithDelay' };
			% set channel, and marker tag character (to postpend to
			% markerbase string to properly identify Marker property name)
			if C == L
				tagchar = 'L';
				obj.Channel = L;
			else
				tagchar = 'R';
				obj.Channel = R;
			end
			% assign values to tags in Stimulus object
			for n = 1:length(markerbase)
				obj.(markerbase{n}) = M.([markerbase{n} tagchar]);
			end

		end	% END buildStimulusFromMarker
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
end	% End of classdef






