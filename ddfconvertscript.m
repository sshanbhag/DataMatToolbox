close all; clear all; clear classes;

%------------------------------------------------------------
%% file names
%------------------------------------------------------------
DLLName = 'C:\DataWave\DWShared\nsDWFile.dll';
datapath = 'F:\Work\Data\BatRestrained\853\';
datapath = 'F:\Work\Data\Bat\BatRestrainedData\FRAdata\';
matpath = 'F:\Work\Data\Bat\BatRestrainedData\FRAdata';
ddffilelist = dir(fullfile(datapath, '*.ddf'));
%save('ddflist.mat', 'ddffilelist');
nDDFfiles = length(ddffilelist);
fprintf('Found %d .DDF file in %s\n\n', nDDFfiles, datapath);

%------------------------------------------------------------
%% convert to matfiles
%------------------------------------------------------------
for n = 1:nDDFfiles
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