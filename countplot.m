Cmatrix	=	[	...
	10	0	11	6	0	0	0	0	2	1	;	...
	0	23	0	2	0	0	1	0	1	3	;	...
	1	0	19	3	0	0	0	0	6	1	;	...
	2	3	3	17	0	0	0	1	2	2	;	...
	5	2	0	0	14	1	6	1	0	1	;	...
	0	0	3	0	0	25	2	0	0	0	;	...
	0	7	0	0	2	4	14	0	3	0	;	...
	0	0	2	5	0	0	0	17	4	2	;	...
	1	0	3	5	0	0	0	2	16	3	;	...
	3	2	1	6	0	0	0	3	4	11	;	...
]

% use pcolor to draw patches with color given by Cmatrix
% % note that y-axis values need to be "flipped" in order for plot to match
% % orientation of Cmatrix
h_plot = pcolor(1:10, 10:-1:1, Cmatrix);

%%
% set the grayscale color map
%%

% # of grayscale levels is set by max value of Cmatrix
ngraylevels = max(max(Cmatrix));
% graymin sets the arbitrary lower bound for the grayscales used
% in the plot.  0 corresponds to black, higher values will raise 
% the minimum value brightness
graymin = 3;
graymax = ngraylevels;

% create grayscale colormap matrix
cvector = (graymin:graymax)./ngraylevels;
cmap = [cvector' cvector' cvector'];

% set the colormap
colormap(cmap);

% draw color bar
h_cb = colorbar;
