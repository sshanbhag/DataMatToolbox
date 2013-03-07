%-----------------------------------------------------------------------------
% Probe.m
%-----------------------------------------------------------------------------
% DataMat Toolbox
% Class Definition
%-----------------------------------------------------------------------------
%	Properties:
% 		t						{nclusters, 1} cell array of spike times.
%								for background data, the format will be
%								{nclusters, nwin} to account for possible
%								multiple background data windows.
% 		wforms				{nclusters, 1} cell array or spike waveforms.  
%								each cluster will be a cell array of nspikes, with
% 								each element containing the sampled waveform
% 		cluster				ID values for units in t and wforms
% 		nclusters			# of clusters (units)
%		samprate				sample rate for waveforms
%		timeunits			units for timestamps in t
%		nwindows				# of spike windows (usually 1, but can be multiple)
% 		name					NeuroShare probe name
%-----------------------------------------------------------------------------
% See also: DWdata, Marker
%-----------------------------------------------------------------------------

%-----------------------------------------------------------------------------
%	Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%-----------------------------------------------------------------------------
% Created: 1 June, 2012 (SJS)
%
% Revisions:
%	25 Jan 2013 (SJS): finally working on this, incorporating NeuroShare
%								bits
%	20 Feb 2013 (SJS): added wform property to store waveforms
%	 -	changed Nclusters to nclusters
%	24 Feb 2013 (SJS): added sample rate property
%	6 Mar 2013 (SJS):
% 	 -	added timeunits property
%	 - add nwindows property
% 	 - added getSpiketimesForCluster method
% 	 - added comments/docs
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
% classdef (ConstructOnLoad = true) Probe < handle
classdef Probe < handle
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Properties
	%------------------------------------------------------------------------
	properties
		t
		wforms
		cluster
		nclusters
		samprate
		time_units
		nwindows
		name
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
		function obj = Probe(varargin)
		%---------------------------------------------------------------------	
		% Probe
		% Constructor method
		%---------------------------------------------------------------------
		% Probe()	
		%							
		%							when called with no arguments, returns empty
		%							Probe object
		%---------------------------------------------------------------------
			DW.DataWaveDefaults;		% load defaults
			%--------------------------------------------------------
			% parse input and verify
			%--------------------------------------------------------
			if isempty(varargin)
				return
			end
		end		% END Probe constructor
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		% General Methods
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		function out = findClusters(obj)
		%---------------------------------------------------------------------	
		% find unique cluster id values
		%---------------------------------------------------------------------
			out = unique(obj.cluster);
			if isempty(obj.nclusters)
				obj.nclusters = out;
			end
		end	% END findClusters
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------

		%---------------------------------------------------------------------
		function spiket = getSpiketimesForCluster(obj, clusternum, varargin)
		%---------------------------------------------------------------------
		% spiket is a [1 X ntimestamps] array of spike timestamps
		%---------------------------------------------------------------------
			% make sure clusternum is inbounds!
			if ~any(clusternum == obj.cluster)
				warning('%s: clusternum %d not found!', mfilename, clusternum)
				fprintf('\tvalid clusters:\t');
				fprintf('%d ', obj.cluster);
				fprintf('\n');
				spiket = [];
				return
			end
			
			if isempty(varargin)
				win = 1;
			else
				if (varargin{1} < 1) || (varargin{1} > obj.nwindows)
					warning('%s: window index %d out of range', ...
																		mfilename, varargin{1});
					fprintf('\tvalid window num: 1 - %d\n', obj.nwindows);
					wmat = [];
					return
				else
					win = varargin{1};
				end
			end
			
			spiket =  cell2mat(obj.t{clusternum == obj.cluster, win});
		end	% END getWaveformsForCluster
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------

		%---------------------------------------------------------------------
		function wmat = getWaveformsForCluster(obj, clusternum, varargin)
		%---------------------------------------------------------------------
		% wmat is a [nwavepoints X ntimestamps] array of spike waveforms
		%---------------------------------------------------------------------
			% make sure clusternum is inbounds!
			if ~any(clusternum == obj.cluster)
				warning('%s: clusternum %d not found!', mfilename, clusternum)
				fprintf('\tvalid clusters:\t');
				fprintf('%d ', obj.cluster);
				fprintf('\n');
				wmat = [];
				return
			end
			if isempty(varargin)
				win = 1;
			else
				if (varargin{1} < 1) || (varargin{1} > obj.nwindows)
					warning('%s: window index %d out of range', ...
																		mfilename, varargin{1});
					fprintf('\tvalid window num: 1 - %d\n', obj.nwindows);
					wmat = [];
					return
				else
					win = varargin{1};
				end
			end
			wmat =  cell2mat(obj.wforms{clusternum == obj.cluster, win})';
		end	% END getWaveformsForCluster
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

