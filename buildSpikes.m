function [Spikes, varargout] = buildSpikes(StimList, varargin)
%------------------------------------------------------------------------
% [Spikes, VarInfo, UnitList] = buildSpikes(StimList, varargin)
%------------------------------------------------------------------------
% DataMat toolbox
%------------------------------------------------------------------------
% build up Spikes cell array to hold spike time information for a 
% given list of Stimulus struct instances
%------------------------------------------------------------------------
% Input Arguments:
%	StimList		array of Stimulus struct objects
% 	
% 	Optional:
% 		UnitList		list of units to include in output Spikes cell array
%
% Output Arguments:
%	Spikes		{Nunits, Nstimuli, Nattenvals} cell array of spiketimes
% 					Spike timestamps are in units of milliseconds (relative to
% 					sweep onset timestamp)
%
%	VarInfo		information about the dimensions of Spikes
% 
% 	UnitList		list of units - really only useful if UnitList is not provided
% 					as an input argument.  In that case UnitList is computed as
% 						[Nunits, tmp] = size(StimList(1).Spiketimes);
% 						UnitList = 1:UnitList;
% 
%------------------------------------------------------------------------
% See also: loadDWfile
%------------------------------------------------------------------------
%
%------------------------------------------------------------------------
% More information:
%------------------------------------------------------------------------
%
% Within each Stimulus (or, StimList in this case) struct element, there are 
% a few key elements for organizing the spikes.
% 
% {R/L}AttenIndices{} contains a list of indices to Stimulus vectors that correspond
% to attenuation values given in {R/L}AttenVals
% 
% for example, 
% 	
% 	The vector at Stimulus(1).RAttenIndices{1} correspondes to a list of 
% 	indices for Stimulus(1) where Stimulus(1)RAttenVals(1) was the 
% 	attenuation setting.  
% 	
% 	So, in order to get the spikestimes corresponding to Stimulus(1) and the 
% 	attenuation value of Stimulus(1).RAttenVals(1), we would use the following
% 	code:
% 	
% 	spikes = Stimulus(1).Spiketimes(<unit number>, Stimulus(1).RAttenIndices{1})
% 		(*** parentheses notation   ^     and                                   ^)
% 		
% 	Which returns a cell array of spike time stamp vectors:
% 	
% 	>> Stimulus(1).Spiketimes(1, Stimulus(1).RAttenIndices{1})
% 
% 		ans = 
% 
% 		 []    [1x15 double]    [1x14 double]
% 
% 	To obtain individual vectors of spike times, use the following notation
% 	(pulling in the second timestamp vector because the first one is empty):
% 	
% 	spikes = Stimulus(1).Spiketimes{1, Stimulus(1).RAttenIndices{1}(2)}
% 	
% 	>>	spikes = Stimulus(1).Spiketimes{1, Stimulus(1).RAttenIndices{1}(2)}
% 
% 		spikes =
% 
% 		Columns 1 through 10
% 
% 		  2868100     2869800     2870000     2870100     2870300     2870500     2870700     2870800     2871000     2871100
% 
% 		Columns 11 through 15
% 
% 		  2871200     2871300     2871500     2871600     2871700
%----------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 	19 July, 2011 (SJS):
%	-	snipped from plotData script
%
% Revisions:
%	2 August, 2011 (SJS):
%	-	added VarInfo as output
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

%-----------------------------------------------------------
% load defaults
%-----------------------------------------------------------
DataWaveDefaults;

%-----------------------------------------------------------
% check inputs
%-----------------------------------------------------------
if nargin == 0
	error('%s: no Stimulus structs given as input!', mfilename);
end

if ~isstruct(StimList)
	error('%s: StimList is not a struct.', mfilename);
else
	Nstimuli = length(StimList);
end

% check if UnitList was provided as additional argument.
if nargin > 1
	% if so, use that as list of units to include in Spikes
	UnitList = varargin{1};
	Nunits = length(UnitList);
else
	% otherwise, use all spikes available
	[Nunits, tmp] = size(StimList(1).Spiketimes);
	UnitList = 1:Nunits;
end


%----------------------------------------------------------------------
% loop through units
%----------------------------------------------------------------------

% unitindex is used to index the Spikes unit (Spikes{unitindex...}) 
% while unit is used to pull the desired units from UnitList
unitindex = 0;
for unit = UnitList
	% increment unitindex counter
	unitindex = unitindex + 1;
	
	% loop through stimuli (this will thus vary across frequencies or wavfiles
	for stimindex = 1:Nstimuli
		
		% loop through different attenuation values
		for attenindex = 1:length(StimList(stimindex).RAttenVals)
			
			% clear the tmpcell variable to avoid confusing results
			clear tmpcell;
			
			% loop through individual sweeps for this attenuation setting
			for sweepindex = 1:length(StimList(stimindex).RAttenIndices{attenindex})
				% To get the spikes for this sweep, we need to access the index that
				% corresponds to the proper attenuation and stimulus combination.
				% These indices are stored in 
				%	Stimulus(1... nstimuli).RattenIndices{1... # att values}(1...nsweeps)
				% 
				% retrieve this index and store it in spikecol
				spikecol = StimList(stimindex).RAttenIndices{attenindex}(sweepindex);
				
				% assign the spike times to the current tmpcell{} element
				tmpcell{sweepindex} = StimList(stimindex).Spiketimes{unit, spikecol};
				% to work with the spike timestamps, need to
				% (1)	subtract off Spikestart(sweepindex) to get times relative to 
				%		start of sweep
				% (2)	then convert from usec to msec
				if ~isempty(tmpcell{sweepindex})
					tmpcell{sweepindex} = 0.001*(tmpcell{sweepindex} -  StimList(stimindex).Sweepstart(spikecol));
				end
			end
			
			% assign the tmpcell spiketimes to the master Spikes cell array
			% Spikes has three dimensions: 
			%		(1) unit index
			%		(2) stimulus index (varying parameter, e.g., wav file, tone freq)
			%		(3) attenuation value
			Spikes{unitindex, stimindex, attenindex} = tmpcell;			
		end	% end of ATTENINDEX
		stiminfo(stimindex).var = StimList(stimindex).Var;
	end	% end of STIMINDEX	
end	% end of UNITINDEX

if nargout > 1
	varargout{1} = UnitList;
end

if nargout > 2
	
	VarInfo(1).Name = 'unit number';
	VarInfo(1).Value = UnitList;
	
	VarInfo(2).Name = 'stim info';
	VarInfo(2).Value = stiminfo;
	
	
	VarInfo(3).Name = 'attenuation value'
	VarInfo(3).Value = atteninfo;
	
	
	
	varargout{2} = VarInfo;
end
