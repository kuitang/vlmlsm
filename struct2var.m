function struct2var(s)
% See http://stackoverflow.com/questions/3470654/how-can-i-move-variables-into-and-out-of-a-structure-akin-to-load-and-save-in-ma
  cellfun(@(n,v) assignin('caller',n,v),fieldnames(s),struct2cell(s));
end
