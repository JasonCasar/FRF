function [startIndices] = determineRuns(matched_data_mat,startTimesMat)
    %Function takes the starting times and returns the index corresponding to
    %the start of each new run (index of matched_data_mat). 
    %
    %INPUTS:
    %   matched_data_mat: matrix containing concentrations, along with the
    %       coordinates and the times they were acquired at. 
    %   startTimesMat: column vector of times corresponding to which time
    %       each individual data run began at (time bin). 
    %
    %OUTPUTS:
    %   startIndices: vector telling us which INDEX in matched_data_mat
    %       corresponds to the start of each individual data run (time bin).

indices = ones(size(startTimesMat,1),1);
for i=1:size(startTimesMat,1)
    startTime = double(datenum(datetime(startTimesMat(i,:),'InputFormat','dd-MMM-yyyy HH:mm:ss')));
    indices(i) = find(matched_data_mat(:,1) == startTime);
end
startIndices = indices;
end
