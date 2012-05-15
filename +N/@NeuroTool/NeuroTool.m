%% ------------------------------------------------------------------------
%++++++++++++++++++++++++++++++++++++++
%++++++This Help File Needs Help+++++++
%++++++++++++++++++++++++++++++++++++++
%   NeuroTool class created to consolidate functions that analyze output
%     of the mat files produced by loadDWfile.m in DataMat Toolbox.
%     (Which are analysis of dw files output by ... [program that records
%     neuro data])
%
%   Authors:Tony Slagle     - tonyslagle@gmail.com
%   Revision history in contents.m
%
%NeuroTool consolidates analysis and handling of interpreted DWFile data.
%
%  dwf = NeuroTool(<converted_DWfile_name>) will load a .mat file that has
%  been converted by loadDWfile.m and prepare it for analysis. If no input
%  then it will open a dialog box asking for input.
%% ------------------------------------------------------------------------
classdef (ConstructOnLoad = true) NeuroTool < handle & N.LoadDWfileOutput
  properties
    psthbinsize;%size of bins. Default is NeuroTool.default('psthbinsize')
  end
%  properties (GetAccess = public, SetAccess = private)
%    psthbins;   %array holding psth bins or maybe this will be a class
%  end
  properties (Access = private)
    psthbinsexist = 0;
  end
  methods
% not yet implemented:
    %ret = binssave(obj, varargin)%saves bins to following naming convention:
                              %folder = sourcefilename
                              %bins files = sourcefilename_probe_unit 
    %a function to catch .psthbins and return a summary of bins calculated
        %by unit, access to units by (#) and a .psthbins.saveCSV to dump to
        %a file
%% ------------------------------------------------------------------------
    function obj = NeuroTool(varargin)
%NeuroTool constructor - dwf = NeuroTool(<filename>,<psthbinsize>)

%parse input args
      fileName = {};
      binSize = N.defaults('psthbinsize');
      arg = 0;  %which arg we're looking at
      while arg < nargin %so we have a block to jump out of
        %if we come to this line twice, we didn't use all the args
        if arg ~= 0
          error('NeuroTool:badinputs','Input args bad. (see help NeuroTool for proper constructor usage)');
        end
        arg = 1;
%filename
        if isequal(class(varargin{arg}), 'char')
          if exist(varargin{arg},'file')
            fileName = varargin(arg);
          else
            warning('NeuroTool:sillyfname','%s does not exist!',varargin{arg});
          end
          arg = arg + 1; %next arg
          if arg > nargin, break; end %jump out of the loop. all args used
        end
%psthbinsize
        if isequal(class(varargin{arg}), 'double')
          binSize = varargin{arg};
          arg = arg + 1; %next arg
          if arg > nargin, break; end %jump out of the loop. all args used
        end
      end%while
      obj = obj@N.LoadDWfileOutput(fileName{:});
      obj.psthbinsize = binSize;
    end
%% ------------------------------------------------------------------------
    function set.psthbinsize(obj,value)
%sets bin size and warns if outside default.psthbsexpected range
      psthexp = N.defaults('psthbsexpected');
      if (value < psthexp(1)) || (value > psthexp(2))
        warning('NeuroTool:psthbinsize:psthbsexpected',...
          'psthbinsize %f outside default range. Time is in seconds.',value)
      end
      obj.psthbinsize = value;
      obj.calcbins;  %recalcualte bins
    end
%% ------------------------------------------------------------------------
    function disp(obj)
      disp@N.LoadDWfileOutput(obj);
      fprintf('\tpsthbinsize:\n\t\t%g s\n\n',obj.psthbinsize)
    end%disp
  end%methods
end%classdef
