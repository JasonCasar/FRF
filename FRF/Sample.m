function [n, obs_index, Z_obs, St, M] = Sample(Z_t,S,ML,MW)
    %This function is nearly identical to the function used by Cressie et
    %al. in their 2010 FRF paper. 
    %It "samples" the entire spatial field (ie, all N bins), returning only
    %those bins for which we have data (ie, either the bins we actually
    %sampled at, or the bins where we haven't removed data from).
    %Currently this is only used in model testing. If it is going to be
    %translated to actual real-time data acquisition, the code needs to be
    %altered to save a pointer to each bin that has a data point in it
    %instead of using ML and MW to store where we are removing data from.
    %
    %INPUTS:
    %   Z_t: The FULL data set (N data points) at time t.
    %   S: The FULL basis function matrix
    %   ML: number telling you which bin to start REMOVING data at.
    %   MW: number telling you how many bins are to be removed. The bins
    %       will be removed from ML to ML+MW (ie, contiguous). 
    %
    %OUTPUTS:
    %   n: number of data points in the removed data set, Z_obs (less than
    %       or equal to N).
    %   obs_index: The index at which there are observations. IE, if we
    %       removed data from bins 10 to 25, then indices 1 to 9 will have a 1,
    %       10 to 25 will have a zero, and 26 to N will have a 1 again.
    %   Z_obs: Continuous vector of data, but with "removed" regions
    %       removed. The length will be less than or equal to N.
    %   St: The basis function matrix of size n<=N x number basis
    %       functions. Again, all regions where data has been "removed" have
    %       been removed from this matrix.
    %   M: Indicator matrix of size Nxn<=N with a 1 everywhere along the 
    %       diagonal where we have data available (ie, not removed).
        
    
N = length(Z_t);
obs_index=ones([N 1]);

last=ML+MW-1;
if last<=N
    obs_index(ML:last,1)=0;
end

Z_obs=Z_t(obs_index>0);
St = S(obs_index(:,1)>0,:);
n = sum(obs_index(:,1));

%Initialize M so that there are zeros everywhere. We will add 1's (ie, add
%an indication that we have data available at that location) at all places
%that are NOT in our removal region. 
M=zeros([N n]);
%If data is being removed this time run (ie, ML>0), we put 1s everywhere data is
%available.
if ML>0
    for ijk=1:N
        %if we haven't gotten to the removal region yet, put a 1 in that
        %diagonal location
        if ijk<ML
              M(ijk,ijk)=1;
        %if we are past the removal region, put a 1 in that OFF-diagonal
        %location. 
        elseif ijk>last
              M(ijk,ijk-MW)=1;
        end
    end
%Otherwise, if data is NOT being removed this time run, we put a 1 in every
%location of the OFF-diagonal. 
else
    for ijk=(last+1):N
        M(ijk,ijk-MW)=1;
    end 
end

end

