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

fname = 'test.txt';
%-----------------------------------------------------------
% count the number of lines in the file
%-----------------------------------------------------------
Nlines = countTextFileLines(fname);
disp(['Found ' num2str(Nlines) ' in file ' fname]);

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
	data = cell(Nlines - 2, nfields);
	size(data)
	
	dcount = 1;
	for n = 2:Nlines-1
		line_n = fgetl(fp);
		tmp = textscan(line_n, '%s', nfields, 'Delimiter', '\t');
		for m = 1:length(tmp{1})
			data{dcount, m} = tmp{1}(m);
		end
		dcount = dcount+1;

	end
end
	

%a = textscan(fp, '%s', 'Delimiter', '\t');


fclose(fp);

% need to post-process the data a bit to place into suitable vectors




