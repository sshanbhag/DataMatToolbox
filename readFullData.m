%% filename
% FullData6.csv is comma-delimited text file exported from FullData6.xlsx
% spreadsheet
datafile = '/Users/sshanbhag/Work/Data/DataWave/batmat/FullData6.csv';
% # of data lines in FullData6 (from looking at xlsx file)
nlines = 8774; 
% matfile in which data will be saved
matfile = '/Users/sshanbhag/Work/Data/DataWave/batmat/FullData6.mat';

%% open file for reading, text mode
fp = fopen(datafile, 'rt');

%% 1st line is list of 1s and 0s where 1 == number and 0 == string or character
% read in line with fgetl, scan comma-delimited text with csvscan, and 
% then convert from cell array to matrix
header.isnumber = str2num(cell2mat(csvscan(fgetl(fp))));
% get # of total columns
ncolumns = length(header.isnumber);

%% read header text
header.fields = csvscan(fgetl(fp));
% fix minor error for last field (which is empty)
header.fields{end + 1} = header.fields{end};

%% need to read in file line-by-line and parse into individual cell values
FullData = cell(nlines, ncolumns);
for n = 1:nlines;
	tmp = textscan(fgetl(fp), '%s', 'Delimiter', ',', 'MultipleDelimsAsOne', 0);
	FullData(n, :) = tmp{1}';
end

%% close file
fclose(fp);

%% save data
save(matfile, '-MAT', 'FullData', 'header', 'nlines', 'ncolumns', 'datafile');