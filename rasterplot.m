function rasterplot(spiketimes, timeMinMax, axesHandle)
%------------------------------------------------------------------------
% H = rasterplot(spiketimes, timeMinMax, axesHandle)
%------------------------------------------------------------------------
% PlotTools toolbox
%------------------------------------------------------------------------
% 
% Given a cell "vector" (NX1 or 1XN) of spiketimes, e.g.
% 
% 			spiketimes = 
% 				 [1x4 double]
% 				 [1x3 double]
% 							  []
% 				 [1x3 double]
% 				 [1x3 double]
%  
% draw a raster plot of the spiketimes, where each element in the spiketimes
% cell vector corresponds to a vector of spiketimes in milliseconds.
% 
% so, for example, spiketimes{1} = [100 120 200 300]
%	 
% timeMinMax is optional x-axis limits as a [1X2] vector; if not provided
% rasterplot will automatically determine limits
% 
% axes is an optional handle to an axes object - if not provided, 
% rasterplot will create a new one in a new figure
% 
%------------------------------------------------------------------------
% Input Arguments:
% 	spiketimes		cell "vector" (NX1 or 1XN) of vectors containing 
% 						times (NOT interspike intervals!) of spikes
% 	
% 	Optional:
% 		timeMinMax		x-axis limit vector in form [min max]
% 		axesHandle		handle to axes
% 
% Output Arguments:
% 	H	handle to plot
%
%------------------------------------------------------------------------
% See also: 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad J. Shanbhag
%	sharad.shanbhag@einstein.yu.edu
%------------------------------------------------------------------------
% Created: 3 December, 2009 (SJS)
%
% Revisions:
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

TICKASCII = double('|');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some checks on inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check if figure handle was provided at input
if exist('axesHandle', 'var')
	% if so, make sure it is a proper handle, if not, create new figure
	if ishandle(axesHandle)
		H = gca(axesHandle);
		
	elseif isempty(axesHandle)
		H = gca;
		
	else
		warning([mfilename ': invalid input axes handle, creating new handle']);
		H = gca;
	end
else
	% otherwise, create axes, save handle to return as output of function
	H = gca;
end

cla

% check if timeMaxMin was specified
if ~exist('timeMinMax', 'var')
	timeMinMax = [0 0];
	maxtimeSearchFlag = 1;
else
	maxtimeSearchFlag = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% draw plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get size of spiketimes
nReps = length(spiketimes);

% loop through reps
for r = 1:nReps
	% x locations for ticks == spike times
	xlocs = spiketimes{r};
	% ylocations are set by rep index (r)
	ylocs = r*ones(size(spiketimes{r}));
	% need a row vector of ticks due to a peculiarity of the text() function
	% in Matlab 
	tickchars = char(TICKASCII * ones(length(spiketimes{r}), 1));
	
	% draw the ticks, return a vector of handles
	h = text(xlocs, ylocs, tickchars);
	% use the handles vector to set color
	set(h, 'Color', [0 0 1]);

	% get the time limit if it wasn't specified by the user
	if maxtimeSearchFlag
		if max(spiketimes{r}) > timeMinMax(2)
			timeMinMax(2) = max(spiketimes{r});
		end
	end
end

% set x limit to manual, set limit
xlim('manual')
xlim(timeMinMax);
% set ylimit to manual, set limit
ylim('manual');
ylim([0 nReps+1]);


	
