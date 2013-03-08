close all; clear all; clear classes;

%------------------------------------------------------------
%% file names
%------------------------------------------------------------
DLLName = 'C:\DataWave\DWShared\nsDWFile.dll';
datapath = 'G:\';
matpath = 'Z:\Marie\batmat';
% ddffilelist = dir(fullfile(datapath, '*.ddf'));
% save('ddflist.mat', 'ddffilelist');

ddffilelist.name = '827_12-18-2012--2954_syllables block (30, 40)_Sorted.ddf';
FORCEFLAG = 1;
%------------------------------------------------------------
%% convert to matfiles
%------------------------------------------------------------
for n = 1:length(ddffilelist)
	ddffile = fullfile(datapath, ddffilelist(n).name);
	[tmppath, tmpfile, tmpext] = fileparts(ddffile);
	matfile = fullfile(matpath, [tmpfile '.mat']);
	if ~exist(matfile, 'file') || FORCEFLAG
		fprintf('Converting:\n');
		fprintf('\t%s\n', ddffile);
		fprintf('\t\t--TO--\n');
		fprintf('\t%s\n\n', matfile);
 		D = DW.convertDDF2MAT(ddffile, 'MATFILE', matfile, 'EVENT', 'SEGMENT', 'NEURAL');
		clear D;
	end
end