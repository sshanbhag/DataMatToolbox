close all

nrows = 4;
ncols = 3;


horizgap = 0.05;
vertgap = 0.05;

plotwidth = (1 - ((ncols+1) * horizgap)) / ncols;
plotheight = (1 - ((nrows+1) * vertgap)) / nrows;

for c = 1:ncols
	xpos(c) = horizgap + ((c-1) * (plotwidth + horizgap));
end

for r = 1:nrows
	ypos(r) = 1 - r*plotheight - r*vertgap;
end

for r = 1:nrows
	for c = 1:ncols
		position{r, c} = [xpos(c) ypos(r) plotwidth plotheight];
	end
end

figure(1)

plotindex = 0;
for r = 1:nrows
	for c = 1:ncols
		plotindex = plotindex + 1;
		
		subplot('Position', position{r, c});
		
		xlim([0 1]);
		ylim([0 1]);
		text(0.05, 0.5, sprintf('r=%d c=%d indx=%d', r, c, plotindex));
		
	end
end

set(gcf, 'Name', 'Subplot Test: single plot')
set(gcf, 'Position', [71   520   560   420]);


horizgap = 0.05;
vertgap = 0.055;
plotgap = 0.0125;

plotwidth = (1 - ((ncols+1) * horizgap)) / ncols;
plotheight = (1 - ((nrows+2) * vertgap)) / (2*nrows);

for r = 1:nrows
	ypos1(r) = 1 - r*plotheight - (r-1)*plotheight - r*vertgap - (r-1)*plotgap;
	ypos2(r) = ypos1(r) - plotheight - plotgap;
	for c = 1:ncols
		xpos(c) = horizgap + ((c-1) * (plotwidth + horizgap));
		pos1{r, c} = [xpos(c) ypos1(r) plotwidth plotheight];
		pos2{r, c} = [xpos(c) ypos2(r) plotwidth plotheight];
	end
end

figure(2)

plotindex = 0;
for r = 1:nrows
	for c = 1:ncols
		plotindex = plotindex + 1;
		
		subplot('Position', pos1{r, c});
		xlim([0 1]);
		ylim([0 1]);
		tmp{1} = sprintf('r=%d c=%d indx=%d', r, c, plotindex);
		tmp{2} = sprintf('x=%.3f  y=%.3f', pos1{r, c}(1), pos1{r, c}(2));
		text(0.05, 0.75, tmp{1});
		text(0.05, 0.25, tmp{2});
		set(gca, 'xtick', []);
		
		if c ~= 1
			set(gca, 'ytick', []);
		end
			
		
		subplot('Position', pos2{r, c});
		xlim([0 1]);
		ylim([0 1]);
		tmp{1} = sprintf('r=%d c=%d', r, c);
		tmp{2} = sprintf('x=%.3f  y=%.3f', pos2{r, c}(1), pos2{r, c}(2));
		text(0.05, 0.75, tmp{1});
		text(0.05, 0.25, tmp{2});
		if r ~= nrows
			set(gca, 'xticklabel', []);
		end
		if c ~= 1
			set(gca, 'ytick', []);
		end
		
	end
end


set(gcf, 'Name', 'Subplot Test: Double (2 row) plot')
set(gcf, 'Position', [640   520   560   420]);


