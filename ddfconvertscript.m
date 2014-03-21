close all; clear all; clear classes;

%------------------------------------------------------------
%% file names
%------------------------------------------------------------
% name of DLL
DLLName = 'C:\DataWave\DWShared\nsDWFile.dll';

% data path is path to input data directory
%datapath = 'F:\Work\Data\BatRestrained\853\';
datapath = 'F:\Work\Data\Bat\BatRestrainedData\BatSorted\866_Sorted\';
% matpath is output path for mat files
matpath = 'F:\Work\Data\Bat\BatRestrainedData\NewMatFiles';
if ~exist(matpath, 'dir')
	error('%s: output dir %s not found', mfilename, matpath)
end
% look for ddf files in data directory
ddffilelist = dir(fullfile(datapath, '*.ddf'));
%save('ddflist.mat', 'ddffilelist');
nDDFfiles = length(ddffilelist);
fprintf('Found %d .DDF file(s) in %s\n\n', nDDFfiles, datapath);

%------------------------------------------------------------
%% convert to matfiles
%------------------------------------------------------------
%for n = 1:nDDFfiles
for n = 7:10
	ddffile = fullfile(datapath, ddffilelist(n).name);
	[tmppath, tmpfile, tmpext] = fileparts(ddffile);
	matfile = fullfile(matpath, [tmpfile '.mat']);
	if ~exist(matfile, 'file')
		fprintf('Converting:\n');
		fprintf('\t%s\n', ddffile);
		fprintf('\t\t--TO--\n');
		fprintf('\t%s\n\n', matfile);
 		D = DW.convertDDF2MAT(ddffile, 'MATFILE', matfile, 'EVENT', 'SEGMENT', 'NEURAL');
		clear D;
	else
		fprintf('File %s already converted....\n\n', ddffile);
	end
	pause
end