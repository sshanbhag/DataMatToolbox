clear all
close all
clear classes

%------------------------------------------------------------
%% file names
%------------------------------------------------------------
Dpath = '/Users/sshanbhag/Work/Data/DataWave/batmat/convertedDDF';
dobjpath = '/Users/sshanbhag/Work/Data/DataWave/batmat/dobj';
unitdatafile = '/Users/sshanbhag/Work/Data/DataWave/batmat/UnitData.mat';
dobjfilesfile = '/Users/sshanbhag/Work/Data/DataWave/batmat/UniqueDobjFiles.mat';

%% load unit data and get file names
load(unitdatafile);

%% find filenames
if ~exist(dobjfilesfile, 'file')
	filenamecol = finddatacolumn('Filename', unitheader.fields);
	% avoid repeats
	dobjfiles = findUniqueUnsortedText(UnitData(:, filenamecol));
else
	disp('Loading dobj files from mat file to save time')
	load(dobjfilesfile);
end

%% process filenames
nfiles = length(dobjfiles)
nfiles = 10;
FORCEFLAG = 1;

for n = 1:nfiles
	% trap situations where filename has " in it
	if dobjfiles{n}(1) == '"'
		dobjfiles{n} = dobjfiles{n}(2:end);
	end
	if dobjfiles{n}(end) == '"'
		dobjfiles{n} = dobjfiles{n}(1:(end-1));
	end
	% parse filename
	[tmppath, tmpname, tmpext] = fileparts(dobjfiles{n});
	% remove Dobj
	tmpname = tmpname(6:end);
	% get unit number
	unitstr = '';
	c = length(tmpname);
	while c >= 1
		if tmpname(c) ~= '_'
			c = c - 1;
		else
			unitstr = tmpname((c+1):end);
			break
		end
	end
	if isempty(unitstr)
		error('Unitstr ERROR!!!!')
	end
	% store unitnumber
	unitnum(n) = str2num(unitstr);
	% store basename
	basename{n} = tmpname(1:c-1);
	fprintf('basename(%d) = %s\n', n, basename{n})
end

% build matDfile list and make sure files exist
for n = 1:nfiles
	matDfile{n} = fullfile(Dpath, [basename{n} '.mat']);
	if exist(matDfile{n}, 'file')
		fprintf('matDfile(%d) = %s\t|\tExists = %d\n', n, basename{n}, exist(matDfile{n}, 'file'))
	else
		error('matDfile(%d) = %s NOT FOUND!', n, matDfile{n})
	end
end


%------------------------------------------------------------
%% create Data Object
%------------------------------------------------------------
clear DW.RateData
probenum = 1;

for n = 1:nfiles
	clear d
	if ~exist(fullfile(dobjpath, dobjfiles{n}), 'file') || FORCEFLAG
		% load matDfile to get D
		fprintf('%s: loading .mat file %s\n', mfilename, matDfile{n});
		load(matDfile{n});
	
		%convert to obj
		% construct d from RateData class
		d = DW.RateData(D, matDfile{n});
		% find vars
		d.findVarsAndAtten;
		% plot unit waveforms (overlaid)
		% d.plotUnitWaveforms('probe', probenum, 'unit', unitnum(n));	

		% save Data Object
		fprintf('Writing d object to %s...\n', fullfile(dobjpath, dobjfiles{n}));
	 	save( fullfile(dobjpath, dobjfiles{n}), 'd', '-MAT');
		
		close all
	end
end
