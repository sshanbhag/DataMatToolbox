function [H, plotopts] = rasterpsthmatrix(Spikes, Nrows, Ncols, plotopts)
%------------------------------------------------------------------------
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
%		ticksymbol		symbol for raster ticks (default is '|')
% 		axesHandle		handle to axes
% 
% Output Arguments:
% 	H	handle to plot
%
%------------------------------------------------------------------------
% See also: rasterplot, psth
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: %	7 July, 2011 (SJS)
%
% Revisions:
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------


% plot a raster and psth for each unit
if ~exist('plotopts', 'var')
	plotopts = struct( ...
		'time_limits',		[0 1000]		, ...
		'horizgap',			0.05			, ...
		'vertgap',			0.055			, ...
		'plotgap',			0.0125		...
	);
end

plotwidth = (1 - ((Ncols+1) * plotopts.horizgap)) / Ncols;
plotheight = (1 - ((Nrows+2) * plotopts.vertgap)) / (2*Nrows);

pos1 = cell(Nrows, Ncols);
pos2 = cell(Nrows, Ncols);

for r = 1:Nrows
	ypos1(r) = 1 - r*plotheight - (r-1)*plotheight - r*plotopts.vertgap - (r-1)*plotopts.plotgap;
	ypos2(r) = ypos1(r) - plotheight - plotopts.plotgap;
	for c = 1:Ncols
		xpos(c) = plotopts.horizgap + ((c-1) * (plotwidth + plotopts.horizgap));
		pos1{r, c} = [xpos(c) ypos1(r) plotwidth plotheight];
		pos2{r, c} = [xpos(c) ypos2(r) plotwidth plotheight];
	end
end

handles1 = cell(Nrows, Ncols);
handles2 = cell(Nrows, Ncols);


for row = 1:Nrows
	for col = 1:Ncols
		% select subplot location for rasters (pos1)
		subplot('Position', pos1{row, col});
		% plot raster
		% store the axes handle returned by rasterplot in the handles2 cell array
		handles1{row, col} = rasterplot(Spikes{row, col}, plotopts.time_limits, '.', 14);
		% turn off xtick labels, and turn off yticks
		set(gca, 'XTickLabel', []);
		set(gca, 'ytick', []);
	
		if isfield(plotopts, 'idlabel')
 			if (col == 1) && (row == 1)
				idstr = plotopts.idlabel;
			else
				idstr = '';
			end
		end
		
		if isfield(plotopts, 'columnlabels')
 			if (row == 1)
				colstr = plotopts.columnlabels{col};
			else
				colstr = '';
			end
		end
			
		if isfield(plotopts, 'rowlabels')
			if col == 1
				rowstr = plotopts.rowlabels{row};
			else 
				rowstr = '';
			end
		end
		
		titlestr = [idstr ' ' colstr ];
		title(titlestr, 'Interpreter', 'none')
		ylabel(rowstr, 'Interpreter', 'none')
	
		% select subplot location for psth (pos2)
		subplot('Position', pos2{row, col});
		% build psth from spike data and plot using bar() function
		[histvals, bins] = psth(Spikes{row, col}, 5, 1000);
		% store the axes handle returned by bar in the handles2 cell array
		handles2{row, col} = bar(bins, histvals);
		% update time limits to match raster
		xlim(plotopts.time_limits)
		% turn off x tick labels for all but the bottom row and
		% turn off y tick labels for all but the left column
		if row ~= Nrows
			set(gca, 'XTickLabel', []);
		end
		if col ~= 1
			set(gca, 'ytick', []);
		end
		set(gca, 'Box', 'off')

		% label the x-axis 'msec' on the lower left psth plot
		if (col == 1) && (row == Nrows)
			xlabel('msec')
		end
		% label the lower right plot axes with the input data file 
		if (col == Ncols) && (row == Nrows)  && isfield(plotopts, 'filelabel')
			xlabel(plotopts.filelabel, 'Interpreter', 'none');
		end

	end
end

H{1} = handles1;
H{2} = handles2;

plotopts.plotwidth = plotwidth;
plotopts.plotheight = plotheight;
plotopts.pos1 = pos1;
plotopts.pos2 = pos2;

		

