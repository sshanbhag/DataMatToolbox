
%------------------------------------------------------------
%% file names
%------------------------------------------------------------
DLLName = 'C:\DataWave\DWShared\nsDWFile.dll';
datapath = 'E:\Bat - SingleCh Restrained\batmat';
dobjpath = 'E:\Bat - SingleCh Restrained\batmat\matfiles';
matfilelist = dir(fullfile(datapath, '*.mat'));

%------------------------------------------------------------
%% create Data Object
%------------------------------------------------------------
for n = 1%:length(matfilelist)
	
	probenum = 1;
	unitnum = 255;
	
	matfile = fullfile(datapath, matfilelist(n).name);
	[tmppath, tmpfile, tmpext] = fileparts(matfile);
	dobjfile = fullfile(dobjpath, [tmpfile '.mat']);

	% load matfile to get D 
	fprintf('%s: loading .mat file %s\n', mfilename, fullfile(datapath, matfile));
	load(matfile);
	
	d = DW.RateData(D, fullfile(datapath, matfile));
	% plot unit waveforms (overlaid)
	d.plotUnitWaveforms('probe', probenum, 'unit', unitnum);
	

	%------------------------------------------------------------
	%------------------------------------------------------------
	%% save Data Object
	%------------------------------------------------------------
	%------------------------------------------------------------
	% write the data object.  for debugging use, the D struct 
	% is also saved, but this will ultimately be redundant...
	save( fullfile(datapath, ['Dobj_' matfile]), 'd', 'probenum', 'unitnum');

	
	
% 	if ~exist(dobjfile, 'file')
% 		fprintf('Converting:\n');
% 		fprintf('\t%s\n', matfile);
% 		fprintf('\t\t--TO--\n');
% 		fprintf('\t%s\n\n', dobjfile);
%  		
% 		clear D;
% 	end
end
