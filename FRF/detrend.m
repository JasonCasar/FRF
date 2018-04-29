function [trend,detrended_data] = detrend(central_coord,data,fullSet,missingStart,missingEnd)
%detrend current run of data. fullSet is a boolean input indicating whether
%or not the regression should be performed on the full data set or only a
%portion of it (ie, the portion where the data are available). 
%If fullSet is set to false, then missingStart and missingEnd are required.
%
%INPUTS:
%   central_coord: Vector of coordinates representing bin locations. Used
%       for calculating regression line.
%   data: The binned concentration measurements
%   fullSet: Boolean. Are all data available or not?
%   missingStart: Only necessary paramter if fullSet==false. Tells us which
%       bins to start removing data at.
%   missingEnd: Again, only necessary if fullSet==false. Tells us which
%       bins to stop removing data at. 
%
%OUTPUTS:
%   trend: Vector containing a value for each bin. This is the trend line
%       calculated from the available data from this time bin.
%   detrended_data: data-trend

data_temp = data;

%Make the data that are "removed" NaN
if fullSet == false ;
    data_temp(missingStart:missingEnd) = NaN;
end
%Perform a linear regression on all AVAILABLE data to generate a trend line
    X0 = [ones(length(central_coord),1) central_coord];
    betahat = regress(data_temp',X0);
    trend=(X0*betahat)';
%Save that trend and subtract it from the data. 
    detrended_data = (data-trend)';
end
