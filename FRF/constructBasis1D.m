function [S,eta,resolution] = constructBasis1D(central_coord,type,r,resolution,signal)
    %Need to add functionality for other basis function types
    if type == 'Gaussian' 
        len = r*(1-2^resolution)/(-1); %resolution of 3 would be len = r + 2r + 4r
        while len>length(central_coord)
           warning('Resolution is being decreased so that the number of basis functions does not exceed the number of bins. Consider decreasing resolution and increasing rank')
           resolution = resolution-1;
           len = r*(1-2^resolution)/(-1);
        end
        if resolution == 0
           error('0 Resolution. Program will crash now. r must be at most smaller than numBins.')
        end
        S = zeros([length(central_coord) len]);
        pointer = 1;
        for res=1:resolution
            numBasisFunc = r*(2^(res-1));
            minCoord = central_coord(1);
            maxCoord = central_coord(length(central_coord));
            range=maxCoord-minCoord;
            spacing=ceil(range/numBasisFunc);
            remainder=mod(range,spacing);
            centers=zeros([1 numBasisFunc]);
            x=minCoord+(remainder/2);
            i=1;
            while i<=numBasisFunc
                centers(1,i)=x;
                x=x+spacing;
                i=i+1;
            end
            sigma=spacing/2; %This is pretty arbitrary
            distsq=zeros([length(central_coord) numBasisFunc]);
            for i=1:numBasisFunc
                distsq(:,i)=(central_coord-centers(i)).^2;
            end
            S_res = exp(-0.5*distsq./(sigma^2));
            S(:,pointer:(pointer+numBasisFunc-1)) = S_res;
            pointer = pointer + numBasisFunc;
        end
        eta = S\signal;
    else
        warning('Currently Only Type Gaussian is Supported At This Point')
    end
end
