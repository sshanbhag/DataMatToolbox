%-----------------------------------------------------------------------------
% Now, we have found the unique STIMULUS TYPES in the stimuli presented
% for each channel.  However, it is still not known if there are differences 
% within the stimulus types; e.g., attenuation values may differ between 
% stimuli, frequency bandwith might vary, etc.  This next chunk of code will
% attempt to sort this out
%-----------------------------------------------------------------------------


% build a matrix of stimulus parameters


% loop through unique indices
for n = 1:Nunique
	% 
	