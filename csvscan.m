function [out, nread] = csvscan(in, varargin)
%------------------------------------------------------------------------
% [out, nread] = csvscan(in, multdelim)
%------------------------------------------------------------------------
% 
% scans comma-delimited text data IN and returns OUT cell vector
% setting multdelim to 1 will treat multiple delimiters to 1
% 
%------------------------------------------------------------------------
% Input Arguments:
% 	Input		input info
% 
% Output Arguments:
% 	Output	output info
%
%------------------------------------------------------------------------
% See also: 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 2013 (SJS)
%
% Revisions:
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

if ~isempty(varargin)
	multdelim = varargin{1};
else
	multdelim = 0;
end

tmp = textscan(	in, '%s', ...
						'Delimiter', ',', ...
						'MultipleDelimsAsOne', multdelim);
out = tmp{1};
nread = length(out);
clear tmp
