%function 
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

time_limits = [0 1000];

Natten = length(StimList(1).RAttenVals);

nrows = Natten;
ncols = Nstimuli;

horizgap = 0.05;
vertgap = 0.055;
plotgap = 0.0125;

plotwidth = (1 - ((ncols+1) * horizgap)) / ncols;
plotheight = (1 - ((nrows+2) * vertgap)) / (2*nrows);

pos1 = cell(nrows, ncols);
pos2 = cell(nrows, ncols);

for r = 1:nrows
	ypos1(r) = 1 - r*plotheight - (r-1)*plotheight - r*vertgap - (r-1)*plotgap;
	ypos2(r) = ypos1(r) - plotheight - plotgap;
	for c = 1:ncols
		xpos(c) = horizgap + ((c-1) * (plotwidth + horizgap));
		pos1{r, c} = [xpos(c) ypos1(r) plotwidth plotheight];
		pos2{r, c} = [xpos(c) ypos2(r) plotwidth plotheight];
	end
end

handles1 = cell(D.Info.Nunits, nrows, ncols);
handles2 = cell(D.Info.Nunits, nrows, ncols);

% for unit = 1:D.Info.Nunits
for unit = 1:1
	% for each unit, create a new page
	H = figure(unit);

	for atten = 1:Natten
		for var = 1:Nstimuli
			% select subplot location for rasters (pos1)
			subplot('Position', pos1{atten, var});
			% plot raster
			% store the axes handle returned by rasterplot in the handles2 cell array
			handles1{unit, atten, var} = rasterplot(Spikes{unit, var, atten}, time_limits);
			% turn off xtick labels, and turn off yticks
			set(gca, 'XTickLabel', []);
			set(gca, 'ytick', []);
			
			if (var == 1) & (atten == 1)
				unitstr = sprintf('Unit %d: ', unit);
			else
				unitstr = '';
			end

			if (atten == 1)
				if strcmp(StimList(var).Type, 'TONE')
					stimstr = sprintf('Stim = %.2f', StimList(var).Var(1).values(1));
				elseif strcmp(StimList(var).Type, 'WAVFILE')
					tmpstr = textscan(StimList(var).Var(2).values{1}, '%s', 'Delimiter', '\\');
					stimstr = ['Stim = ' tmpstr{1}{end}];
				else
					stimstr = ['Stim = ' num2str(StimList(var).Var(1).values(1))];
				end
			else
				stimstr = '';
			end
			
			if var == 1
				attenstr = sprintf('%d dB', StimList(var).RAttenVals(atten));
			else 
				attenstr = '';
			end
				
			titlestr = [unitstr ' ' stimstr ];
			title(titlestr, 'Interpreter', 'none')
			ylabel(attenstr, 'Interpreter', 'none')
	
			% select subplot location for psth (pos2)
			subplot('Position', pos2{atten, var});
			% build psth from spike data and plot using bar() function
			[histvals, bins] = psth(Spikes{unit, var, atten}, 5, 1000);
			% store the axes handle returned by bar in the handles2 cell array
			handles2{unit, atten, var} = bar(bins, histvals);
			% update time limits to match raster
			xlim(time_limits)
			% turn off x tick labels for all but the bottom row and
			% turn off y tick labels for all but the left column
			if atten ~= Natten
				set(gca, 'XTickLabel', []);
			end
			if var ~= 1
				set(gca, 'ytick', []);
			end
			
			% label the x-axis 'msec' on the lower left psth plot
			if (var == 1) && (atten == Natten)
				xlabel('msec')
			end
			% label the lower right plot axes with the input data file 
			if (var == Nstimuli) && (atten == Natten)
				xlabel(D.Info.file, 'Interpreter', 'none');
			end
				
		end
	end
end

		

