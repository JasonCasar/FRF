function [W,matched_data_mat,central_coord,count] = matchBin(numBins,time_nmea,xcoord,ycoord,time_mcpc,aveconc,startTimesMat)
    %Given data from the gps and the mcpc, match all data points that occur at
    %the same time, and bin these data points into numBins bins. This is
    %1D. 

    
    %INPUTS: 
    %   numBins - numbins - This will be N, the number of evenly spaced locations 
    %        you predict at.
    %   time_nmea - array of the times each x,y pair were recorded at
    %   xcoord - array of the x coordinates from nmea
    %   ycoord - array of the y coordinates from nmea
    %   time_mcpc - array of the times each conc data were recorded at
    %   aveconc - array of the concentration values from mcpc
    %   startTimesMat - Column vector indicating when each time starts
    %        (includes training sets). 
    %
    %OUTPUTS:
    %   W -  a CELL containing number of runs vectors (number of different start times) 
    %       by numBins vectors. Each vector contains the coordinates
    %       pointing to the location in matched_data_mat where you can find
    %       the measuredment that occured in that bin.
    %   matched_data_mat - The measurement in each space-time bin
    %   central_coord - The central coordinate representing each bin
    %   count - number of data points per bin

% Find all instances of time_mcpc that are in time_nmea (the gps has a
% faster acquisition rate, so every member of time_mcpc should appear in
% time_nmea). 
[members,Loc] = ismember(time_mcpc,time_nmea);
time_mcpc = datenum(time_mcpc);
matched_data_mat = [];
for i=1:length(members)
    if members(i)==1;
        % combine these corresponding measurements into a matrix that
        % includes a single time, xcoordination, ycoordinate and
        % concentration.
        matched_data = [time_mcpc(i) xcoord(Loc(i)) ycoord(Loc(i)) aveconc(i)];
        matched_data_mat = vertcat(matched_data_mat,matched_data);
    end
end

%call determineRuns helper function to translate startTimes into start
%indices in matched_data_mat (ie, what indices does that time correspond
%to?)
startIndices = determineRuns(matched_data_mat,startTimesMat);
%Remove all data from matched_data_mat that does not occur within the first
%start time and the last (end) time (ie, we haven't started driving yet, or
%we've finished sampling and are driving back but haven't turned off the
%instruments yet). 
matched_data_mat = matched_data_mat(startIndices(1):startIndices(length(startIndices)),:);
%Because we've removed those indices we now need to make all startIndices
%relative to the first start index (startIndices(1). 
startIndices = startIndices - startIndices(1)+1;

%Save all coordinates that occur in the first time run
firstCoords = matched_data_mat(startIndices(1):startIndices(2),2);
%minCoord is the minimum latitude coordinate. It's possible that we started
%sampling at the minimum or maximum coordinate, so we need to determine
%which one we did. 
if firstCoords(1)<firstCoords(length(firstCoords));
   minCoord = firstCoords(1);
   maxCoord = firstCoords(length(firstCoords));
else
    minCoord = firstCoords(length(firstCoords));
    maxCoord = firstCoords(1);
end
%MinCoord and maxCoord are the extremes of our sampling field.
%width of each bin is the length of this field divided by the number of
%bins we want. 
binWidth = abs(maxCoord-minCoord)/numBins;

%Initiate an empty cell that is number of times by number of spatial bins.
%Each of those TxN bins will contain a vector of arbitrary size that
%contains all indices of the data in matched_data_mat that fall within that
%particular space-time bin pair. 
W = cell((length(startIndices)-1),numBins);
%Initialize an empty vector that will contain the central location of each
%bin.
central_coord = zeros(1,numBins);
%intialize empty matrix that will contain the number of data points that
%fall within that space-time bin pair. 
count = zeros(length(startIndices)-1,numBins);
for i=1:(length(startIndices)-1)
    %Determine all coordinates that are in this time run
    thisRunCoords = matched_data_mat(startIndices(i):startIndices(i+1),2);
    j=1;
    %Loop through all coordinates between minCoord and maxCoord, 
    %incrementing by binWidth.
    for coordinate=minCoord:binWidth:maxCoord
        central_coord(j) = coordinate+(binWidth/2);
        %Find the index of all data that fall within that time period and
        %spatial bin. Place those indices in W so that you can reference
        %them later. 
        W{i,j} = find(abs(thisRunCoords-coordinate)<(binWidth/2))+startIndices(i)-1;
        count(i,j) = length(W{i,j});
        j=j+1;
    end
end
%if you want to know what W looks like, uncomment the line below. 
%celldisp(W)
end


