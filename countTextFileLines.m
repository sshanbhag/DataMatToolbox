function [nLines, nCols] = countTextFileLines(filename)
%------------------------------------------------------------------------
% nLines = countTextFileLines(filename)
%------------------------------------------------------------------------
% 
% if <filename> is a text file, counts number of lines in it
%
% returns 0 if file is binary or on error
% 
%------------------------------------------------------------------------
% Input Arguments:
% 	filename		file name
% 
% Output Arguments:
% 	nLines		# of lines in text file.
%					0 on error
%
%------------------------------------------------------------------------
% See also: fopen, fgetl, fclose
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neoucom.edu
%------------------------------------------------------------------------
% Created: 24 January, 2011 (SJS)
%
% Revisions:
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------


%----------------------------------------------------------------
% check that file exists; warn, and return 0 if file not found
%----------------------------------------------------------------
if ~exist(filename, 'file')
	warning('MATLAB:file', 'mfilename: file %s not found', mfilename, filename);
	nLines = 0;
	return
end

%----------------------------------------------------------------
% open file, throw error if not opened successfully
%----------------------------------------------------------------
[fp, msg] = fopen(filename, 'rt');
if fp == -1
	error('%s: error opening file <%s> (%s)', mfilename, filename, msg);
end

%----------------------------------------------------------------
% get file information
%----------------------------------------------------------------
[FILENAME,PERMISSION,MACHINEFORMAT,ENCODING] = fopen(fp);

%----------------------------------------------------------------
% Count lines by looping through the file
%----------------------------------------------------------------
nLines = 0;
while ~feof(fp)
	tmp = fgetl(fp);
	nLines = nLines+1;
end

if nargout == 2
	%----------------------------------------------------------------
	% count columns for each line
	%----------------------------------------------------------------
	% preallocate nCols to hold column counts per line
	nCols = zeros(nLines, 1);
	% rewind file pointer to beginning of file
	frewind(fp);
	% loop through lines, scan each line of text, save # of elements
	for n = 1:nLines
		tmpln = fgetl(fp);
		tmp = textscan(tmpln, '%s', 'Delimiter', '\t');
		nCols(n) = length(tmp{1});
	end
end
%----------------------------------------------------------------
% close file
%----------------------------------------------------------------
fclose(fp);
