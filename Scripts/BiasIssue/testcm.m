% reset random # gen
RandObj = RandStream.getDefaultStream;
RandObj.reset;

Stimuli = 2:15;
Reps = 5:5:100;
Iterations = 100;

% pre allocate arrays
Hmean = zeros(length(Stimuli), length(Reps));
Hsd = zeros(length(Stimuli), length(Reps));

% loop through nstimuli list
for s = 1:length(Stimuli)
	Nstimuli = Stimuli(s);
	% loop through nreps list
	for r = 1:length(Reps)
		Nreps = Reps(r);
		
		H = zeros(1, Iterations);
		% loop for # of iterations
		for n = 1:Iterations
			% create matrix of counts for confusion matrix
			% confusion matrix has # stimuli along rows, test along column
			%
			% need to generate using monte carlo method (to ensure that 
			% votes along each row add up to consistent value)

			CM = zeros(Nstimuli);
			% loop through rows (responses)
			for cm_s = 1:Nstimuli
				% loop for number of trials
				for cm_r = 1:Nreps
					rbin = randi(Nstimuli, 1);
					CM(cm_s, rbin) = CM(cm_s, rbin) + 1;
				end
			end
			H(n) = cminfo(CM);
		end	% END Iterations Loop
		Hmean(s, r) = mean(H);
		Hsd(s, r) = std(H);
	end	% END Nresp Loop
end	% END Nstimuli Loop

display_matrix(Hmean, '%.2f');

display_matrix(Hsd, '%.2f');


corder = [	0.00	0.00	1.00	;	...
				0.00	1.00	0.00	;	...
				1.00	0.00	0.00	;	...
				0.00	0.25	0.75	;	...
				0.00	0.75	0.25	;	...
				0.25	0.75	0.00	;	...
				0.75	0.25	0.00	;	...
				0.00	0.50	0.50	;	...
				0.50	0.00	0.50	;	...
				1.00	0.00	1.00	;	...
				1.00	1.00	0.00	;	...
				0.00	1.00	1.00	;	...
				0.00	0.00	0.00	;	...
			];

figure(1)
set(gcf, 'DefaultAxesColorOrder', corder);
subplot(211)
plot(Reps, Hmean', '.-')
legend(num2str(Stimuli'))
legend(gca, 'Location', 'Best')
legend boxoff
ylabel('Hmean (bits)')
title('Random Confusion Matrix')

subplot(212)
plot(Reps, Hsd', '.-');
xlabel('number of trials per stimulus')
ylabel('H std. dev. (bits)')
title('colors indicate # of stimuli')