function [W,data_detrend,X0,betahat,central_coord,count] = matchBin(numBins,time_nmea,xcoord,ycoord,time_mcpc,aveconc,startTimesMat)
%[W,data_detrend,X0,betahat,central_coord,count] = matchBin(10,time_nmea,xcoord,ycoord,time_mcpc,aveconc,startTimesMat);

%Given data from the gps and the mcpc, match all data points that occur at
%the same time, and bin these data points into numBins bins. This is 1D
%W will be a cell containing number of runs vectors by numBins vectors.
%this function also calls data_detrend. so the data it returns will be
%detrended.
%Each vector will contain the coordinates pointing to the location in
%data_detrend where you can find the measurement that occured in that bin.
%To find the coordinate of that bin, look in central_coord.

[members,Loc] = ismember(time_mcpc,time_nmea);
time_mcpc = datenum(time_mcpc);
matched_data_mat = [];
for i=1:length(members)
    if members(i)==1;
        matched_data = [time_mcpc(i) xcoord(Loc(i)) ycoord(Loc(i)) aveconc(i)];
        matched_data_mat = vertcat(matched_data_mat,matched_data);
    end
end

%call determineRuns helper function
startIndices = determineRuns(matched_data_mat,startTimesMat);

%detrend. Call detrend helper function
[data_detrend,X0,betahat] = detrend(matched_data_mat,startIndices);
startIndices = startIndices - startIndices(1)+1;

firstCoords = data_detrend(startIndices(1):startIndices(2),2);
if firstCoords(1)<firstCoords(length(firstCoords));
   minCoord = firstCoords(1);
   maxCoord = firstCoords(length(firstCoords));
else
    minCoord = firstCoords(length(firstCoords));
    maxCoord = firstCoords(1);
end
binWidth = abs(maxCoord-minCoord)/numBins;

W = cell((length(startIndices)-1),numBins);
central_coord = zeros(1,numBins);
count = zeros(length(startIndices)-1,numBins);
for i=1:(length(startIndices)-1)
    thisRunCoords = data_detrend(startIndices(i):startIndices(i+1),2);
    j=1;
    for coordinate=minCoord:binWidth:maxCoord
        central_coord(j) = coordinate+(binWidth/2);
        W{i,j} = find(abs(thisRunCoords-coordinate)<(binWidth/2))+startIndices(i)-1;
        count(i,j) = length(W{i,j});
        j=j+1;
    end
end
%celldisp(full_data)
end


