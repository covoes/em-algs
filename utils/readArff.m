function [data]=readArff(file)
%READARFF reads the data in an ARFF file and put in a matlab matrix, only can be used in unix (egrep command needed)
data = [];

cmd = sprintf('egrep -v "(^(@|%%)|^\\s*$)"  %s > temp', file);
[s,r] = system(cmd);

data = csvread('temp');
