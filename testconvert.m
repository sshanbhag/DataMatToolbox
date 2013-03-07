close all; clear all; clear classes;

%------------------------------------------------------------
%% file names
%------------------------------------------------------------
DLLName = 'C:\DataWave\DWShared\nsDWFile.dll';
datapath = 'G:\';
matpath = 'Z:\Marie\batmat';
ddffilelist = dir(fullfile(datapath, '*.ddf'));
save('ddflist.mat', 'ddffilelist');

%------------------------------------------------------------
%% convert to matfiles
%------------------------------------------------------------
for n = 1:length(ddffilelist)
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
	end
end