function Spikes = build_simdata(varargin)
%------------------------------------------------------------------------------
% Spikes = build_simdata(varargin)
%---------------------------------------------------------------------------
% 
% builds fake spiketrain datasets
% 
%---------------------------------------------------------------------------
% Input Arguments:
% 
% Output Arguments:
%
%---------------------------------------------------------------------------
% See also: 
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%---------------------------------------------------------------------------
% Created: 4 April, 2013 (SJS)
%
% Revisions:
%---------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
% check inputs
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
% set default values
set_defaults;

if ~isempty(varargin)
	argn = 1;
	while argn <= length(varargin)
		switch upper(varargin{argn})
			case 'DEFAULT'
				break;
			case 'CATEGORIES'
				Ncategories = varargin{argn+1};
				argn = argn + 2;
			case 'REPS'
				Nreps = varargin{argn+1};
				argn = argn + 2;
			case 'NSPIKES'
				Nspikes = varargin{argn+1};
				argn = argn + 2;
			case 'TOFFSET'
				Toffset = varargin{argn+1};
				argn = argn + 2;
			case 'RANDOM_ISI'
				RandISI = 1;
				MeanISI = varargin{argn+1};
				argn = argn + 2;
			case 'RANDOM_NSPIKES'
				RandNspikes = 1;
				MeanNspikes = varargin{argn+1};
				argn = argn + 2;
			otherwise
				error('%s: unknown option %s', mfilename, varargin{argn});
		end	% END switch(varargin)
	end	% END while argn
else
	fprintf('%s: using default values\n', mfilename);
end	% END if

if length(Nspikes) == 1
	Nspikes = Nspikes * ones(Ncategories, 1);
end

if RandNspikes && (length(MeanNspikes) == 1)
	MeanNspikes = MeanNspikes * ones(Ncategories, 1);
end


%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
% create array
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
Spikes = cell(Ncategories, Nreps);
if RandISI
	for c = 1:Ncategories
		for r = 1:Nreps
			Spikes{c, r} = poisson_spiketrain(Nspikes(c), MeanISI);
		end
	end
else
	for c = 1:Ncategories
		base_train = Toffset(c):(1/Nspikes(c)):1;
		for r = 1:Nreps
			Spikes{c, r} = base_train;
		end
	end
end
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
end
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------


%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
function set_defaults
	evalin('caller', 'Ncategories = 3;');
	evalin('caller', 'Nreps = 5;');
	evalin('caller', 'Nspikes = [5 5 5];');
	evalin('caller', 'Toffset = [.01 .02 .03];');
	evalin('caller', 'RandISI = 0;');
	evalin('caller', 'MeanISI = 1;');
	evalin('caller', 'RandNspikes = 0;');
	evalin('caller', 'MeanNspikes = 5;');
end
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
