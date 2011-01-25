% function Output = function_template(Input)
%------------------------------------------------------------------------
% Output = function_template(Input)
%------------------------------------------------------------------------
% 
% Description
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
% sshanbhag@neoucom.edu
%------------------------------------------------------------------------
% Created: 24 January, 2011 (SJS)
%
% Revisions:
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------



fp = fopen('test.txt', 'r');

header.line1 = fgetl(fp);
header.line2 = fgetl(fp);

l1fields = textscan(header.line1, '%s', 'Delimiter', '\t');

% scan the second header line for strings (tab delimited)
% these strings indicate the data in successive rows
l2fields = textscan(header.line2, '%s', 'Delimiter', '\t')
nfields = length(l2fields{1});

if ~nfields
	warning('DWFILE:empty-data', ...
				'%s: no fields found in header line 2 of data file', mfilename);
else
	n = 1;
	STOP_SCAN = 0;
	
	while ~feof(fp) && ~STOP_SCAN
		line_n = fgetl(fp);
		tmp = textscan(line_n, '%s', nfields, 'Delimiter', '\t');
		if feof(fp)
			STOP_SCAN = 1;
		else
			n = n+1;
		end
	end
end
	

%a = textscan(fp, '%s', 'Delimiter', '\t');


fclose(fp);




