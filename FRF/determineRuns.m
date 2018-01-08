function [startIndices] = determineRuns(matched_data_mat,startTimesMat)
%Function takes the starting times and returns the index corresponding to the start of each new run. 
indices = ones(size(startTimesMat,1),1);
for i=1:size(startTimesMat,1)
    startTime = double(datenum(datetime(startTimesMat(i,:),'InputFormat','dd-MMM-yyyy HH:mm:ss')));
    indices(i) = find(matched_data_mat(:,1) == startTime);
end
startIndices = indices;
end
