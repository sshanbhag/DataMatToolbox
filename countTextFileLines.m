function nLines = countTextFileLines(filename)
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


if ~exist(filename, 'file')
	warning('MATLAB:file', 'mfilename: file %s not found', mfilename, filename);
	nLines = 0;
	return
end


[fp, msg] = fopen(filename, 'rt');

if fp == -1
	error('%s: error opening file <%s> (%s)', mfilename, filename, msg);
end

[FILENAME,PERMISSION,MACHINEFORMAT,ENCODING] = fopen(fp);

nLines = 0;

while ~feof(fp)
	tmp = fgetl(fp);
	nLines = nLines+1;
end

fclose(fp);


