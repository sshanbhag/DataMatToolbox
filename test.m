nrows = 3;
ncols = 4;


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

posindex = 1;
for r = 1:nrows
	for c = 1:ncols
		position{posindex} = [xpos(c) ypos(r) plotwidth plotheight];
		posindex = posindex + 1;
	end
end









plotindex = 0;
for r = 1:nrows
	for c = 1:ncols
		plotindex = plotindex + 1;
		
		subplot('Position', position{plotindex});
		
		title(num2str([r c plotindex]));
	end
end

		