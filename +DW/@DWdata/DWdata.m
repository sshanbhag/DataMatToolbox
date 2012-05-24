%-----------------------------------------------------------------------------
% DWdata.m
%-----------------------------------------------------------------------------
% DataMat Toolbox
% Class Definition
%-----------------------------------------------------------------------------
%	DWdata facilitates access to data output by loadDWfile.m
%
%  obj = DWdata('xxx_y_zz_ttt_n.mat') loads the following:
%     obj.fpath;   %path to file
%     obj.fname;   %name of file
%     obj.fext;    %file extension
%     obj.Background; %really this should be a class too...
%     obj.D;          %really this should be a class too...
%     obj.Stimulus;   %really this should be a class too...
%     obj.fullfname;  %full file name
%-----------------------------------------------------------------------------
% See also: loadDWfile (function)
%-----------------------------------------------------------------------------

%-----------------------------------------------------------------------------
% initial coding and design:
%	Tony Slagle
%	tonyslagle@gmail.com
% Continuing development: 
%	Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%-----------------------------------------------------------------------------
% Created: January, 2012 (TS)
%
% Revisions:
%	17 May, 2012 
%		updated documentation
%		renamed file and class from LoadDWfileData.m to DWdata.m
%-----------------------------------------------------------------------------
% TO DO:
%-----------------------------------------------------------------------------

classdef (ConstructOnLoad = true) DWdata < handle
	
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Define protected properties
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	properties (SetAccess = protected)
		fpath;			%path to file
		fname;			%name of file
		fext;			%file extension
		Background;	%really this should be a class too...
		D;				%really this should be a class too...
		Stimulus;		%really this should be a class too...
	end
	
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Define other properties
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	properties (Dependent = true)
		fullfname;  %full file name
	end
  
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Define methods
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	methods
	
		
	%% ------------------------------------------------------------------------
	function obj = DWdata(varargin)
	%---------------------------------------------------------------------	
	%	DWdata(<fileName>) 
	%	Constructor method
	%	opens file called fileName (char) or opens
	%  a dialog box to get a filename if the fileName provided does not exist.
	%---------------------------------------------------------------------	

		%parse input and verify
		obj.fname = '';
		if nargin > 1
			error('DWdata:toomanyinputs','too many inputs! Try DWdata(filename)');
		elseif nargin == 1
			if exist(varargin{1},'file') == 2
				[obj.fpath, obj.fname, obj.fext] = fileparts(varargin{1});
				if isempty(obj.fpath)
					obj.fpath = pwd;
				end
			else
				warning('DWdata:sillyfname','%s does not exist!',varargin{1});
			end
		end

		%if still don't have a file name, get one
		if isequal(obj.fname,'')
			[fileName, obj.fpath] = uigetfile('*.mat','Open .mat file outupt from loadDWfile');
			if fileName == 0
				error('DWdata:nochoose','no filename chosen')
			end
			[~, obj.fname, obj.fext] = fileparts(fileName);
		end

		%Try loading the file and see if crashes
		try
			loaded = load(obj.fullfname, 'D', 'Stimulus', 'Background');
			obj.D = loaded.D;
			obj.Stimulus = loaded.Stimulus;
			obj.Background = loaded.Background;
		catch
			error('%s: could not find file %s', mfilename, obj.fullfname);
		end
	end		%DWdata
	
%% ------------------------------------------------------------------------
	function disp(obj)
	%------------------------------------------------------
	% displays information about loaded data
	%------------------------------------------------------
		fprintf(1, '\t%s - Loaded from file:\n',class(obj));
		fprintf(1,'\t\t%s\n\n',obj.fullfname);
		fprintf(1,'Contains data:\n');
		fprintf('\tStimulus:\n\t\t%ix%i %s\n',size(obj.Stimulus,1),size(obj.Stimulus,2),class(obj.Stimulus))
		fprintf('\tBackground:\n\t\t%ix%i %s\n',size(obj.Background,1),size(obj.Background,2),class(obj.Background))
		fprintf('\tD:\n\t\t%ix%i %s\n',size(obj.D,1),size(obj.D,2),class(obj.D))
	end	%disp
	
%% ------------------------------------------------------------------------
	function ret = get.fullfname(obj)
		ret = fullfile(obj.fpath, [obj.fname obj.fext]);
	end	%get.fullfname
	
	
	% End of methods
	end
end