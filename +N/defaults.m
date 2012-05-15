function ret = defaults(varargin)
%N.defaults allows saving NeuroTool presets like PSTH bin size, etc
% N.defaults returns the NeuroTool defaults stored in NTdefaults.mat:
%     .psthbinsize    default bin size to group spikes into
%     .psthbsexpected expected bin size range
%
% ret = NeuroTool.defaults(<varname>) returns the value of varname
%
% NeuroTool.defaults(<varname>, <value>) sets varname to value
%
% See also N

  % Added 11 Jan 2012 by SJS to work with MAC platform
  if any(strcmpi(computer, {'GLNX86', 'GLNXA64', 'MACI64'}))
	  load('/+N/NTdefaults.mat')
  else
	  load('+N\NTdefaults.mat')
  end
  
  switch nargin
    case 0
      ret = NTdefaults;
    case 1
      ret = NTdefaults.(varargin{1});
    case 2
      NTdefaults.(varargin{1}) = varargin{2};
      ret = NTdefaults;
      % Added 11 Jan 2012 by SJS to work with MAC platform
      if any(strcmpi(computer, {'GLNX86', 'GLNXA64', 'MACI64'}))
        save('+N/NTdefaults.mat','NTdefaults')
      else
        save('+N\NTdefaults.mat','NTdefaults')
      end
    otherwise
      error('NTdefaults:toomanyargs','Too many inputs, see NeuroTool.defaults help for usage.');
  end
end