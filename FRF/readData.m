function [xcoord,ycoord,time_nmea,aveconc,time_mcpc] = readData(nmea_file,mcpc_file,startTimesMat)
    %returns an array for the x and y coordinates, an array for concentation 
    %and time arrays for nmea and mcpc.
    %
    %INPUTS: 
    %      IMPORTANT: CHANGE dayof, monthof and yearof below
    %      nmea_file - A text file containing coordinates of measurements
    %      mcpc_file - The raw text output of the mcpc 
    %      startTimesMat - Column vector indicating when each time starts
    %        (includes training sets) 
    %
    %OUTPUTS:
    %      xcoord - array of the x coordinates from nmea
    %      ycoord - array of the y coordinates from nmea
    %      time_nmea - array of the times each x,y pair were recorded at
    %      time_mcpc - array of the times each conc data were recorded at
    %      aveconc - array of the concentration values from mcpc
    

%%%%% THESE MUST BE CHANGED %%%%
dayof = 8;
monthof = 12;
yearof = 2017;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileID = fopen(nmea_file);
A = textscan(fileID,'%s');
fclose(fileID);
A=A{1};

xcoord = strjoin(A(1:4:length(A)));
xcoord = strsplit(xcoord,{'(',', (',','});
xcoord = xcoord(2:(length(xcoord)-1));
xcoord = cellfun(@str2num,xcoord);
xcoord = xcoord';

ycoord = strjoin(A(2:4:length(A)));
ycoord = strsplit(ycoord,{')',' '});
ycoord = ycoord(1:(length(ycoord)-1));
ycoord = cellfun(@str2num,ycoord);
ycoord = ycoord';

time_nmea_cell = strjoin(A(4:4:length(A)));
time_nmea_cell = strsplit(time_nmea_cell,' ');
time_nmea = datetime(time_nmea_cell,'InputFormat','HH:mm:ss');
time_nmea = time_nmea';
dayofvec = dayof*ones(length(time_nmea),1);
monthofvec = monthof*ones(length(time_nmea),1);
yearofvec = yearof*ones(length(time_nmea),1);
time_nmea.Year = yearofvec;
time_nmea.Month = monthofvec;
time_nmea.Day = dayofvec;

fileID = fopen(mcpc_file);
A = textscan(fileID,'%s');
fclose(fileID);
A=A{1};
A=A(38:length(A));

aveconc = A(3:21:length(A));
aveconc = cellfun(@str2num,aveconc);

time_mcpc_cell = A(2:21:length(A));
time_mcpc = datetime(time_mcpc_cell,'InputFormat','HH:mm:ss');
dayofvec = dayof*ones(length(time_mcpc),1);
monthofvec = monthof*ones(length(time_mcpc),1);
yearofvec = yearof*ones(length(time_mcpc),1);
time_mcpc.Year = yearofvec;
time_mcpc.Month = monthofvec;
time_mcpc.Day = dayofvec;

%The gps recrods greenwhich meantime, not local time. correct this by
%subtracting 8 hours
timeOffset = hours(8);
time_nmea = time_nmea - timeOffset;

%downsample the nmea data by averaging every point that occurs in the same
%second (our text file only has second resolution even though the gps actually measures multiple times a second)
[unique_time_nmea,~, ic] = unique(time_nmea,'stable');
aveXcoord = accumarray(ic,xcoord,[],@mean);
aveYcoord = accumarray(ic,ycoord,[],@mean);
xcoord = aveXcoord;
ycoord = aveYcoord;
time_nmea = unique_time_nmea;

end
