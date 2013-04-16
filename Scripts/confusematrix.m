function varargout = confusematrix(SpkTrains, Cost, zVal)

[nStimuli, nReps] = size(SpkTrains);

% initialize Ddistance array
Ddistance= zeros(nStimuli,nStimuli);

%------------------------------------------------------------------------
% compute metric, build confusion matrix
%------------------------------------------------------------------------
% loop for the number of stimuli
for C1 = 1:nStimuli
	%loop for number of stimulus repetition trials;
	for T1 = 1:nReps
		% allocate Dspike
		Dspike = zeros(nStimuli, nReps);
		% for each pair of stimuli, compute spike distance metric for all trials
		for C2 = 1:nStimuli
			for T2=1:nReps
				% call the spkd() function with the 2 spike trains and cost
				Dspike(C2,T2) = (spkd(SpkTrains{C1,T1}, SpkTrains{C2,T2}, Cost))^zVal;
			end	% END T2 loop
		end	% END C2 loop

		% set the distance from this train to itself to 0
		Dspike(C1,T1) = 0;   
		% trap Inf values, set to 0
		[a,b] = find(Dspike==Inf);
		Dspike(a,b) = 0;

		% compute average distance
		Dsum = sum(Dspike,2);
		Dave = zeros(size(Dsum));
		for C3 = 1:nStimuli,
			if C3 == C1
				Dave(C3) = (Dsum(C3)/(nReps-1))^(1/zVal);
			else
				Dave(C3) = (Dsum(C3)/nReps)^(1/zVal);
			end	% END C3==C1 if
		end	% END C3 loop
		% trap Inf values and set to 0
		[a, b] = find(Dave == Inf);
		Dave(a, b) = 0;
		% find indices of minimum avg distance
		[c ,d] = min(Dave);
		% increment the "vote" for the minimum 
		Ddistance(C1,d) = Ddistance(C1,d) + 1;
	end	% END T1 loop
end	% END C1 loop

varargout{1} = Ddistance;