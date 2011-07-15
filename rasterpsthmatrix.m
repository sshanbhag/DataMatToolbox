function [H, plotopts] = rasterpsthmatrix(Spikes, plotopts)
%------------------------------------------------------------------------
% [H, plotopts] = rasterpsthmatrix(Spikes, plotopts)
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
% 			timelimits: [0 1000]
% 			psth_binwidth: 5
% 			raster_tickmarker: '.'
% 			raster_ticksize: 12
% 			horizgap: 0.0500
% 			vertgap: 0.0550
% 			plotgap: 0.0125
% 			filelabel: '768_4_1q_Bat_1_output.txt'
% 			columnlabels: {7x1 cell}
% 			rowlabels: {4x1 cell}
% 			idlabel: 'Unit 11'
% 
% Output Arguments:
% 	H				handle to plot
%	plotopts		plot options structure, with handles updated
% 		timelimits: [0 1000]
% 		psth_binwidth: 5
% 		raster_tickmarker: '.'
% 		raster_ticksize: 12
% 		horizgap: 0.0500
% 		vertgap: 0.0550
% 		plotgap: 0.0125
% 		filelabel: '768_4_1q_Bat_1_output.txt'
% 		columnlabels: {7x1 cell}
% 		rowlabels: {4x1 cell}
% 		idlabel: 'Unit 11'
% 		plotwidth: 0.0857
% 		plotheight: 0.0837
% 		pos1: {4x7 cell}
% 		pos2: {4x7 cell}
% 
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
		'timelimits',				[0 1000]			, ...
		'psth_binwidth',			5					, ...
		'raster_tickmarker',		'|'				, ...
		'raster_ticksize',		10					, ...
		'horizgap',					0.05				, ...
		'vertgap',					0.055				, ...
		'plotgap',					0.0125			...
	);
end

%----------------------------------------------------------------------------
% compute plot widths and plot heights
%----------------------------------------------------------------------------
% logic:
%-----------
% any figure window is measured in dimensionless units, with x-axis (horizontal)
% having range of [0 1] and y-axis (vertical) [0 1].  The upper-left corner of
% the figure window has coordinates (0, 1) and lower-left corner is (1, 0)
% 
% Therefore to have all plots the same width and height, and to have them evenly
% spaced with horizontal gaps of size plotopts.horizgap and vertical gaps of
% plotopts.vertgap, compute the width as follows:
% 
% 	plotwidth = ([window width] - [# of columns + 1] * [horizgap]) / [# of columns]
% 	
% and the height as follows (noting that the # of plots in the vertical
% dimension is 2 * Nrows, since there are 2 plots (raster & psth) per row):
% 
% 	plotheight = ([window height] - [# of rows + 2] * [vertgap]) / 2*[# of rows]
% 	 
%----------------------------------------------------------------------------
plotwidth = (1 - ((Ncols+1) * plotopts.horizgap)) / Ncols;
plotheight = (1 - ((Nrows+2) * plotopts.vertgap)) / (2*Nrows);

%----------------------------------------------------------------------------
% compute positions for rasters (pos1) and psths (pos2)
%----------------------------------------------------------------------------
% logic:
%-----------
%----------------------------------------------------------------------------
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


%----------------------------------------------------------------------------
% In order to scale all the psths to same min max y values, need to
% pre-compute the psths and keep track of min and max values.  do this here
%----------------------------------------------------------------------------
% data will be stored in psth struct with histvals and bins cell matrices
%		histvals{Nrows, Ncols} will hold histogram data
%		bins{Nrows, Ncols} will hold the respective x-axis bins for the histvals
%		maxval(Nrows, Ncols) will hold the max values of histvals for each element
%----------------------------------------------------------------------------

% pre-allocate them here.
psthdata.histvals = cell(Nrows, Ncols);
psthdata.bins = cell(Nrows, Ncols);
psthdata.maxval = zeros(Nrows, Ncols);

% loop through rows and cols of Spikes and plot data
for row = 1:Nrows
	for col = 1:Ncols
		% compute psth

		% build psth from spike data and plot using bar() function
		[tmpvals, tmpbins] = psth(	Spikes{row, col}, ...
											plotopts.psth_binwidth, ...
											plotopts.timelimits(2));
		psthdata.histvals{row, col} = tmpvals;
		psthdata.bins{row, col} = tmpbins;
		psthdata.maxval(row, col) = max(tmpvals);
	end
end

% find the overall maximum value in psthdata.maxval matrix and build ylimit
% vector
psthdata.ylimits = [0 max(max(psthdata.maxval))];


%----------------------------------------------------------------------------
% Now, plot the data
%----------------------------------------------------------------------------

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
		handles1{row, col} = rasterplot(	Spikes{row, col}, ...
													plotopts.timelimits, ...
													plotopts.raster_tickmarker, ...
													plotopts.raster_ticksize);
												
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

		% store the axes handle returned by bar in the handles2 cell array
		handles2{row, col} = bar(psthdata.bins{row, col}, psthdata.histvals{row, col});
		
		% update time limits to match raster
		xlim(plotopts.timelimits)
		% set ylimits to overall value in psthdata
		ylim(psthdata.ylimits);
		
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

	end		% END OF col LOOP
end		% END OF row LOOP

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
	plotopts.psthlimits = psthdata.ylimits;
end

