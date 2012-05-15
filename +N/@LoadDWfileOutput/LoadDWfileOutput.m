classdef (ConstructOnLoad = true) LoadDWfileOutput < handle
%LoadDWfileOutput facilitates access to data output by loadDWfile.m
%  obj = LoadDWfileOutput('xxx_y_zz_ttt_n.mat') loads the following:
%     obj.fpath;   %path to file
%     obj.fname;   %name of file
%     obj.fext;    %file extension
%     obj.Background; %really this should be a class too...
%     obj.D;          %really this should be a class too...
%     obj.Stimulus;   %really this should be a class too...
%     obj.fullfname;  %full file name
  properties (SetAccess = protected)
    fpath;   %path to file
    fname;   %name of file
    fext;    %file extension
    Background; %really this should be a class too...
    D;          %really this should be a class too...
    Stimulus;   %really this should be a class too...
  end
  properties (Dependent = true)
    fullfname;  %full file name
  end
  methods
%% ------------------------------------------------------------------------
    function obj = LoadDWfileOutput(varargin)
%  LoadDWfileOutput(<fileName>) opens file called fileName (char) or opens
%  a dialog box to get a filename if the fileName provided does not exist.

%parse input and verify
      obj.fname = '';
      if nargin > 1
        error('LoadDWfileOutput:toomanyinputs','too many inputs! Try LoadDWfileOutput(filename)');
      elseif nargin == 1
        if exist(varargin{1},'file') == 2
          [obj.fpath, obj.fname, obj.fext] = fileparts(varargin{1});
            if isempty(obj.fpath)
              obj.fpath = pwd;
            end
        else
          warning('LoadDWfileOutput:sillyfname','%s does not exist!',varargin{1});
        end
      end
%if still don't have a file name, get one
      if isequal(obj.fname,'')
        [fileName, obj.fpath] = uigetfile('*.mat','Open .mat file outupt from loadDWfile');
        if fileName == 0
          error('LoadDWfileOutput:nochoose','no filename chosen')
        end
        [~, obj.fname, obj.fext] = fileparts(fileName);
      end
%Try loading the file and see if crashes
      loaded = load(obj.fullfname, 'D', 'Stimulus', 'Background');
      obj.D = loaded.D;
      obj.Stimulus = loaded.Stimulus;
      obj.Background = loaded.Background;
    end%LoadDWfileOutput
%% ------------------------------------------------------------------------
    function disp(obj)
      fprintf(1, '\t%s - Loaded from file:\n',class(obj));
      fprintf(1,'\t\t%s\n\n',obj.fullfname);
      fprintf(1,'Contains data:\n');
      fprintf('\tStimulus:\n\t\t%ix%i %s\n',size(obj.Stimulus,1),size(obj.Stimulus,2),class(obj.Stimulus))
      fprintf('\tBackground:\n\t\t%ix%i %s\n',size(obj.Background,1),size(obj.Background,2),class(obj.Background))
      fprintf('\tD:\n\t\t%ix%i %s\n',size(obj.D,1),size(obj.D,2),class(obj.D))
    end%disp
%% ------------------------------------------------------------------------
      function ret = get.fullfname(obj)
         ret = fullfile(obj.fpath, [obj.fname obj.fext]);
      end%get.fullfname
   end
end