function [S,eta] = constructBasis1D(coords,type,r,resolution,signal)
   %Need to add functionality for other basis function types and multiple
   %resolutions. 
    if type == 'Gaussian' 
        if resolution == 1
            range=coords(length(coords))-coords(1);
            spacing=ceil(range/r);
            remainder=mod(range,spacing);
            centers=zeros([1 r]);
            x=coords(1)+(remainder/2);
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
            distsq=zeros([length(coords),r]);
            for i=1:r
                %I don't need this second for loop. I can just subtract the
                %whole coordinate vector (not 1:range!!! change that) from
                %centers(i) and then .^2
                for j=1:range %Going to need to redo this because my coordinates don't start at 1. 
                    distsq(j,i)=(j-centers(i))^2;
                end
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
