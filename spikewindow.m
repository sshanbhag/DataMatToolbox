% create fake spike times (milliseconds)
s = [100 200 210 220 230 300 500 1000 1100 1900 2001 2020  4002 4050 4060 4070];

% find spikes between 0 and 300 ms

% first, use between function to find the values in s that are in range
valid_s = between(s, 100, 300)

% then, use the logical array returned by between to get the actual values
s_window = s(valid_s)


% now break up spikes into 1000 ms "sweeps"

tmin = 0:1000:4000;
tmax = tmin + 1000;

% loop through windows
for w = 1:length(tmin)
	spikes{w} = s(between(s, tmin(w), tmax(w)))	
end


% alternative method
for w = 1:length(tmin)
	altspikes{w} = s((s >= tmin(w)) & (s < tmax(w)))
end

% OR  use windowspikes function
