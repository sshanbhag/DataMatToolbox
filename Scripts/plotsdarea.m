function varargout = plotsdarea(xdata, ydata, yerr, varargin)
%------------------------------------------------------------------------------
% fillhandle = plotsdarea(xdata, ydata, yerr, varargin)
%---------------------------------------------------------------------------
% 
% This function will fill a region around the ydata using the standard deviation
% in yerr using the Matlab fill command.
% 
% fillhandle is the returned handle to the filled region in the plot.
%
% 
%---------------------------------------------------------------------------
% Input Arguments:
% xpoints= The horizontal data points (ie frequencies). Note length(Upper)
%         must equal Length(lower)and must equal length(xpoints)!
% upper = the upper curve values (data can be less than lower)
% lower = the lower curve values (data can be more than upper)
% color = the color of the filled area 
% edge  = the color around the edge of the filled area
% add   = a flag to add to the current plot or make a new one.
% transparency is a value ranging from 1 for opaque to 0 for invisible for
% the filled color only.
% 
% Output Arguments:
%
%---------------------------------------------------------------------------
% See also: 
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
% Based on jbfill by John A. Bockstage, November 2006 (Matlab Central)
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%
%
%USAGE: [
%This function will fill a region with a color between the two vectors provided
%using the Matlab fill command.
%
%fillhandle is the returned handle to the filled region in the plot.
%xpoints= The horizontal data points (ie frequencies). Note length(Upper)
%         must equal Length(lower)and must equal length(xpoints)!
%upper = the upper curve values (data can be less than lower)
%lower = the lower curve values (data can be more than upper)
%color = the color of the filled area 
%edge  = the color around the edge of the filled area
%transparency is a value ranging from 1 for opaque to 0 for invisible for
%the filled color only.
%
%John A. Bockstege November 2006;
%Example:
%     a=rand(1,20);%Vector of random data
%     b=a+2*rand(1,20);%2nd vector of data points;
%     x=1:20;%horizontal vector
%     [ph,msg]=jbfill(x,a,b,rand(1,3),rand(1,3),0,rand(1,1))
%     grid on
%     legend('Datr')
%---------------------------------------------------------------------------
% Created: 16 April, 2013 (SJS)
%
% Revisions:
%------------------------------------------------------------------------------
	
	%---------------------------------------------------------------------------
	%---------------------------------------------------------------------------
	% check inputs
	%---------------------------------------------------------------------------
	%---------------------------------------------------------------------------
	% set default values
	set_defaults;

	if ~isempty(varargin)
		argn = 1;
		while argn <= length(varargin)
			switch upper(varargin{argn})
				case 'DEFAULT'
					break;
				case 'TRANSPARENCY'
					transparency = varargin{argn+1};
					argn = argn + 2;
				case 'EDGE'
					edge = varargin{argn+1};
					argn = argn + 2;
				case 'COLOR'
					color = varargin{argn+1};
					argn = argn + 2;
				case 'LINESTYLE'
					linestyle = varargin{argn+1};
					argn = argn + 2;
				otherwise
					error('%s: unknown option %s', mfilename, varargin{argn});
			end	% END switch(varargin)
		end	% END while argn
	end	% END if

	% ensure all data are in row vectors
	if iscolumn(xdata)
		xdata = xdata';
	end
	if iscolumn(ydata)
		ydata = ydata';
	end
	if iscolumn(yerr)
		yerr = yerr';
	end
	
	% check lengths
	nx = length(xdata);
	ny = length(ydata);
	nysd = length(yerr);
	if ~all(nx == [ny nysd])
		msg = 'Error: Must use the same number of points in each vector';
	end

	%---------------------------------------------------------------------------
	%---------------------------------------------------------------------------
	% plot
	%---------------------------------------------------------------------------
	%---------------------------------------------------------------------------
	% upper part of filled error is bounded by ydata + yerr
	yupper = ydata + yerr;
	% lower part is bounded by ydata - yerr
	ylower = ydata - yerr;
	% construct area boundaries
	filled = [yupper, fliplr(ylower)];
	xfilldata = [xdata, fliplr(xdata)];

	% plot the data - hold status is not important, since if it is not set 
	% outside the function, this will overwrite existing plot
	plot(xdata, ydata, linestyle);

	% now, check hold state and enable it if it is not set
	orig_holdstate = ishold;
	if ~orig_holdstate
		hold on
	end

	% create filled area
	fillhandle = fill(xfilldata, filled, color);

	% set edge properties and color
	set(fillhandle, ...	
				'EdgeColor', edge, ...
				'FaceAlpha', transparency, ...
				'EdgeAlpha', transparency	);

	% if orig_holdstate was turned on by the program,
	% (orig_holdstate == 0) turn it off
	if ~orig_holdstate
		hold off
	end
	
	if nargout
		varargout{1} = fillhandle;
	end

end	% END FUNCTION


%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
function set_defaults
	% default is to have a transparency of .5
	evalin('caller', 'transparency = 0.5;');
	% default edge color is black
	evalin('caller', 'edge = ''k'';');	
	% default color is blue
	evalin('caller', 'color = ''b'';');	
	% plot line style
	evalin('caller', 'linestyle = [color ''-''];');
end
