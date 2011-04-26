% function out = teststimulus(duration_ms, Fs, varargin)
%------------------------------------------------------------------------
% Output = function_template(Input)
%------------------------------------------------------------------------
% 
% Description
% 
%------------------------------------------------------------------------
% Input Arguments:
% 	Input		input info
% 
% Output Arguments:
% 	Output	output info
%
%------------------------------------------------------------------------
% See also: 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neoucom.edu
%------------------------------------------------------------------------
% Created: 13 April, 2011 (SJS)
%
% Revisions:
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

Fs = 100000;
duration_ms = 100;

% % compute length of output stimulus
% N = 0.001 * duration_ms * Fs;
% 
% % pre-allocate vector
% out = zeros(N, 1);
% 
% % get the stimulus type
% if ~nargin
% 	type = 'default'
% else
% 	type = lower(varargin{1});
% end


% create canonical "spike"
tlim = 8.65372789;
t = tlim * .001:1/Fs:.001;
K = besselj(0, abs(tlim*t));

[tmp, limindices] = find((t>=-1) & (t<=1));

% spike = sin2array(K(limindices), 1, Fs);
% 
% plot(t(limindices), K(limindices), t(limindices), spike, 'r-.')

plot(t(limindices), K(limindices))

[A, B] = min(abs(0-K))

t(B)
K(B)

% % act depending on stimulus type
% switch type
% 	case 'default'
% 		
% 	otherwise
% 		error('%s: unknown stimulus type %s', mfilename, type);
% end
% 
	