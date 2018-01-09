function [S,eta] = constructBasis1D(central_coord,type,r,resolution,signal)
   %Need to add functionality for other basis function types and multiple
   %resolutions. 
    if type == 'Gaussian' 
        if resolution == 1
            minCoord = central_coord(1);
            maxCoord = central_coord(length(central_coord));
            range=maxCoord-minCoord;
            spacing=ceil(range/r);
            remainder=mod(range,spacing);
            centers=zeros([1 r]);
            x=minCoord+(remainder/2);
            i=1;
            while i<=r
                centers(1,i)=x;
                x=x+spacing;
                i=i+1;
            end
            
            %Adding multiple resolutions would be a simple matter of making
            %more sigmas and then adding a column to the S matrix for each
            %of the new basis functions from the new resolutions. 
            
            sigma=spacing/2; %might change this
            distsq=zeros([length(central_coord),r]);
            for i=1:r
                distsq(:,i)=(central_coord-centers(i)).^2;
            end
            S = exp(-0.5*distsq./(sigma^2));
            eta=S\signal;
            
        else
            print('Currently Only Resolution 1 is Supported At This Point')
        end
    else
        print('Currently Only Type Gaussian is Supported At This Point')

    end
end
