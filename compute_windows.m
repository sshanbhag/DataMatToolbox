function [windows, nwindows] = compute_windows(total_time, window_time)
%-----------------------------------------------------------------------------
% windows = compute_windows(total_time, window_time)
%-----------------------------------------------------------------------------
% 
% computes windows ([start end]) given total time and size of windows
% 
%-----------------------------------------------------------------------------
% Input Arguments:
%
% 	total_time		total time to break into windows
% 	window_time		length of windows
%-----------------------------------------------------------------------------
% Output Arguments:
% 	windows			cell array of 1X2 windows
%
%-----------------------------------------------------------------------------
% See also: 
%-----------------------------------------------------------------------------

%-----------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%-----------------------------------------------------------------------------
% Created: 10 February, 2012 (SJS)
%
% Revisions:
%-----------------------------------------------------------------------------
% TO DO:
%-----------------------------------------------------------------------------

nwindows = round(total_time / window_time);
windows = cell(nwindows, 1);
for b = 1:nwindows
	windows{b} = [b-1 b] * window_time;
end