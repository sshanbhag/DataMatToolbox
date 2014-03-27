function [tag, value] = readTaggedField(S)

% find first ':'
semicol = find(S ==':', 1);

tag = S(1:(semicol - 1));
value = S((semicol + 1):end);



