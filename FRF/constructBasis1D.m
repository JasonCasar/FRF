function [S,eta,resolution] = constructBasis1D(central_coord,type,r,resolution,signal)
    %This function will take in the location of each basis function (of
    %resolution 1), ie central_coord, and return the matrix of basis
    %functions (currently, the only type of basis function supported in
    %gaussian but you can change that to add another type of function if
    %you want to, like a wavelet). It will also take in the signal of data
    %you want to model the initial basis function coefficients off of. 
    %
    %INPUTS:
    %   central_coord: vector telling us where to cener each of the basis
    %       functions of resolution 1.
    %   type: the type of function you want (currently must be gaussian)
    %   r: The number of basis function in the first resolution. if
    %       resolution = 1, this will be the total number of basis functions.
    %   resolution: Resolution 1 will contain r gaussians, each with a
    %       variance of sigma. Resolution 2 will contain 2r gaussians, each
    %       with a variance of sigma/2. Resolution 3 will contain 4r gaussians,
    %       each with a variance of sigma/4. Resolution 4 will contain 8r
    %       gaussians each with a variance of sigma/8 etc. etc. 
    %   signal: The training set that you want to model your basis
    %       coefficients off of. These will become your initial basis function
    %       coefficients. 
    
    
    %OUTPUTS:
    %   S: Nxtotal-number-of-bins matrix of gaussian basis functions 
    %       centered at central_coord
    %   eta: your intial basis function coefficients. Determined from
    %       signal. (ie, S*eta will equal signal).
    %   resolution: just returns what resolution you inputted (don't really
    %       remember why I did this, oh well...)
    
    %Must be of type "Guassian". Capitalization matters. 
    if type == 'Gaussian' 
        len = r*(1-2^resolution)/(-1); %resolution of 3 would be len = r + 2r + 4r
        %This will automatically decrease your resolution until the total
        %number of basis functions is less than or equal to the number of
        while len>length(central_coord)
           warning('Resolution is being decreased so that the number of basis functions does not exceed the number of bins. Consider decreasing resolution and increasing rank')
           resolution = resolution-1;
           len = r*(1-2^resolution)/(-1);
        end
        %If r>numBins, no amount of resolution reduction will make
        %r<numBins, so program will throw the error below. 
        if resolution == 0
           error('0 Resolution. Program will crash now. r must be at most smaller than numBins.')
        end
        %S is Nxtotal number of bins
        S = zeros([length(central_coord) len]);
        pointer = 1;
        %Loop through each resolution and make basis functions for each
        for res=1:resolution
            %Figure out how many basis functions of this resolution there
            %will be
            numBasisFunc = r*(2^(res-1));
            minCoord = central_coord(1);
            maxCoord = central_coord(length(central_coord));
            range=maxCoord-minCoord;
            %Figure out spacing based on numBasisFunc (the higher the
            %resolution, the more basis functions and the closer together
            %they will be). 
            spacing=ceil(range/numBasisFunc);
            remainder=mod(range,spacing);
            centers=zeros([1 numBasisFunc]);
            %Figure our the location of the first basis function
            x=minCoord+(remainder/2);
            i=1;
            %Figure out where all the basis functions are centered and add
            %to centers. 
            while i<=numBasisFunc
                centers(1,i)=x;
                x=x+spacing;
                i=i+1;
            end
            sigma=spacing/2; %This is a pretty arbitrary choice for variance.
            distsq=zeros([length(central_coord) numBasisFunc]);
            %Calculate distance squared of each central coord from all
            %centers, and square that distance. This will be a matrix.
            for i=1:numBasisFunc
                distsq(:,i)=(central_coord-centers(i)).^2;
            end
            %Once you have the squared distances for that central coord
            %from all centers, use it to calculate the gaussian according
            %to the formula for a gaussian function. 
            S_res = exp(-0.5*distsq./(sigma^2));
            %for resolution 1, pointer is zero. then it's 0+r, then
            %it's r+2r, then its 3r+4r etc. 
            S(:,pointer:(pointer+numBasisFunc-1)) = S_res;
            %Once you've gotten all the way through the basis functions in
            %this resolution, increment pointer to next batch of basis
            %functions (ie, in next resolution). 
            pointer = pointer + numBasisFunc;
        end
        %Once you have calculated the entire S matrix, calculate eta, the
        %coefficeints, from S divided by signal. 
        eta = S\signal;
    else
        %You entered a type other than "Gaussian"
        warning('Currently Only Type Gaussian is Supported At This Point')
    end
end
