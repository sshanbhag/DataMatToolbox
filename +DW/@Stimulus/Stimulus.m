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

% class definition
classdef (ConstructOnLoad = true) Stimulus < handle
	%% Properties
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Define protected properties
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	properties
		MarkerList
		Type
		Channel
		Indices
		Nreps
		Var
		Nsweeps
		Sweepstart
		Sweepend
		PreSweep
		PostSweep
		LAttenVals
		LAttenIndices
		RAttenVals
		RAttenIndices
		Spiketimes
		FirstTimestamp
 		Timestamp
		Tagstring
		
		% general properties
		Amplitude
		TimeShift
		RampUp
		RampDown
	end	% end of properties
	
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Define methods
	%------------------------------------------------------------------------
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

			%--------------------------------------------------------
			%parse input and verify
			%--------------------------------------------------------
			if isempty(varargin)
				return
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
		function obj = buildStimulusFromMarkers(obj)
		%---------------------------------------------------------------------	


		end	% END buildStimulusFromMarkers
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






