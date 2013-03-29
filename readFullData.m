% readFullData.m
%
% Script to read in csv file FullData6.xlsx, convert it to a cell matrix, and
% store in FullData6.mat file.  
%
% Additionally, (added 20 Mar 2013) create a UnitInfo matrix of only the 
% unit information, without any statistics, and store in UnitInfo.mat file.
% This can then be used for selecting units/Dobj files for other analyses
% without having to read in all of the FullData matrix values (time consuming!)
%

%% filename
% FullData6.csv is comma-delimited text file exported from FullData6.xlsx
% spreadsheet
datafile = '/Users/sshanbhag/Work/Data/DataWave/batmat/FullData6.csv';
% # of data lines in FullData6 (from looking at xlsx file)
nlines = 8774; 
% matfile in which data will be saved
fullmatfile = '/Users/sshanbhag/Work/Data/DataWave/batmat/FullData6.mat';
unitmatfile = '/Users/sshanbhag/Work/Data/DataWave/batmat/UnitInfo.mat';

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
save(fullmatfile, '-MAT', 'FullData', 'header', 'nlines', 'ncolumns', 'datafile');

%% create UnitData (information-only for FullData)

% find column that has Esigthreshold
ecol = finddatacolumn('Esigthreshold', header.fields);
% unitInfo will hold 1:(ecol - 1) fields
UnitInfo = FullData(:, 1:(ecol - 1));
% unitheader will hold header info
unitheader.fields = header.fields(1:(ecol - 1));
unitheader.isnumber = header.isnumber(1:(ecol - 1))

% save to file
save(unitmatfile, '-MAT', 'UnitInfo', 'unitheader');