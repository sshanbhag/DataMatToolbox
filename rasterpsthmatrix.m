function [H, plotopts] = rasterpsthmatrix(Spikes, Nrows, Ncols, plotopts)
%------------------------------------------------------------------------
% [H, plotopts] = rasterpsthmatrix(Spikes, Nrows, Ncols, plotopts)
%------------------------------------------------------------------------
% PlotTools toolbox
%------------------------------------------------------------------------
% 
% 
%------------------------------------------------------------------------
% Input Arguments:
% 	Spikes		Cell matrix of {Nrows, Ncols}, each element of which contains
% 					a cell array of spike time vectors (NOT interspike intervals!)
% 					in milliseconds 
% 
%	Optional inputs:
%		plotopts		Plot options structure
% 			time_limits: [0 1000]
% 			horizgap: 0.0500
% 			vertgap: 0.0550
% 			plotgap: 0.0125
% 			filelabel: '768_4_1q_Bat_1_output.txt'
% 			columnlabels: {7x1 cell}
% 			rowlabels: {4x1 cell}
% 			idlabel: 'Unit 1'
% 
% Output Arguments:
% 	H				handle to plot
%	plotopts		plot options structure, with handles updated
% 		time_limits: [0 1000]
% 		horizgap: 0.0500
% 		vertgap: 0.0550
% 		plotgap: 0.0125
% 		filelabel: '768_4_1q_Bat_1_output.txt'
% 		columnlabels: {7x1 cell}
% 		rowlabels: {4x1 cell}
% 		idlabel: 'Unit 1'
% 		plotwidth: 0.0857
% 		plotheight: 0.0837
% 		pos1: {4x7 cell}
% 		pos2: {4x7 cell}
%------------------------------------------------------------------------
% See also: rasterplot, psth
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 	7 July, 2011 (SJS)
%
% Revisions:
%	13 July, 2011 (SJS)
% 		-	functionalized script
% 		-	removed Nrows and Ncols as inputs (redundant)
% 		-	updated comments/documentation
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

% get dimensions of Spikes cell matrix
[Nrows, Ncols] = size(Spikes);

% check size
if ~Nrows || ~Ncols
	error('%s: error in size of Spikes cell matrix', mfilename)
end

% initialize plotopts struct if not passed in
if ~exist('plotopts', 'var')
	plotopts = struct( ...
		'time_limits',		[0 1000]		, ...
		'horizgap',			0.05			, ...
		'vertgap',			0.055			, ...
		'plotgap',			0.0125		...
	);
end

% compute plot widths and plot heights
plotwidth = (1 - ((Ncols+1) * plotopts.horizgap)) / Ncols;
plotheight = (1 - ((Nrows+2) * plotopts.vertgap)) / (2*Nrows);

% compute positions for rasters (pos1) and psths (pos2)
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

% initialize cell matrices for storing handles to raster plots (handles1) and 
% psths (handles2)
handles1 = cell(Nrows, Ncols);
handles2 = cell(Nrows, Ncols);

% loop through rows and cols of Spikes and plot data
for row = 1:Nrows
	for col = 1:Ncols

		%-------------------------------------------------------
		% First, plot raster for this row and column
		%-------------------------------------------------------

		% select subplot location for rasters (pos1)
		subplot('Position', pos1{row, col});
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
	
		%-------------------------------------------------------
		% then, plot psth
		%-------------------------------------------------------

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

%-------------------------------------------------------
% set output values
%-------------------------------------------------------

H{1} = handles1;
H{2} = handles2;

if nargout == 2
	plotopts.plotwidth = plotwidth;
	plotopts.plotheight = plotheight;
	plotopts.pos1 = pos1;
	plotopts.pos2 = pos2;
end

