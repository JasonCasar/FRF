function [Y_pred,var_pred,diff,transpose_data,initial_diff] = FixedRankFiltering(nmea_file,mcpc_file,startTimesMat,missing,r,numBins,resolution)
    %FOR STATIC FIXED RANK FILTERING
    %This function performs fixed rank filtering with static parameters
    %Relies on the following helper functions: Sample, constructBasis1D,
    %readData, matchBin, determineRuns.
    %
    %INPUTS:
    %      nmea_file - A text file containing coordinates of measurements.
    %      Needs to be in the same folder as this function.
    %      mcpc_file - The raw text output of the mcpc. Needs to be in the
    %      same folder as this function. 
    %      startTimesMat - Column vector indicating when each time starts
    %        (includes training sets)
    %      missing - number between 0 and 0.99 indicating % of the
    %        field is to be removed for certain time periods. Used for
    %        testing algorithm performance
    %      r - The number of basis functions in the FIRST resolution (if
    %        resolution =1 this is just the number of basis functions)
    %      resolution - number of different resolutions for basis
    %        functions. Each subsequent resolution has 2 times as many
    %        basis functions as the previous one, and a variance for the
    %        half the value of the previous ones.
    %      numbins - This will be N, the number of evenly spaced locations 
    %        you predict at. IMPORTANT: must be chosen so that each bin
    %        contains at least 1 measurement.
    %
    %OUTPUTS:
    %      Y_pred - NxT matrix with predictions at each bin for all T time
    %        periods
    %      var_pred - NxT matrix with associated prediction errors
    %      diff - NxT matrix with differences between predicted and
    %        measured values (used for algorithm testing)
    %      transpose_data - NxT matrix with BINNED meaured data
    %      initial_diff - Nx1 contains differences between two training
    %        sets (used for algorithm testing)
    
    
    %Below is a sample way to run the code for the morning data from Dec
    %8th 2017, with 60 bins and a total of 60 basis functions.
    %%%IMPORTANT:  num basis functions MUST be less than or equal to num bins
    %startTimesMat = ['08-Dec-2017 10:37:41';'08-Dec-2017 10:47:44';'08-Dec-2017 10:57:01';'08-Dec-2017 11:06:56';'08-Dec-2017 11:15:56';'08-Dec-2017 11:25:57';'08-Dec-2017 11:35:00'];
    %nmea_file = 'coordinates.txt'; mcpc_file='MCPC_171208_102434.txt';numBins=60;r=20;resolution=2;missing=0.1;
    %[Y_pred,var_pred,diff,transpose_data] = FixedRankFiltering(nmea_file,mcpc_file,startTimesMat,missing,r,numBins,resolution);
    
    %Below is a sample way to run the code for the afternoon data from Dec
    %8th 2017. 
    %startTimesMat = ['08-Dec-2017 14:01:40';'08-Dec-2017 14:13:39';'08-Dec-2017 14:26:47';'08-Dec-2017 14:36:13';'08-Dec-2017 14:47:31';'08-Dec-2017 14:59:23';'08-Dec-2017 15:12:40'];
    %nmea_file = 'coordinates1.txt'; mcpc_file='MCPC_171208_135234.txt';numBins=60;r=20;resolution=2;missing=0.1;
    %[Y_pred,var_pred,diff,transpose_data] = FixedRankFiltering(nmea_file,mcpc_file,startTimesMat,missing,r,numBins,resolution);    

    %%%%% CODE BEGINS HERE %%%%%%
    
    %Using readData.m, get arrays from the nmea for each x,y pair and when
    %that pair was recorded, and arrays from the mcpc for concentration
    %values and when that concentration was recorded. 
    [xcoord,ycoord,time_nmea,aveconc,time_mcpc] = readData(nmea_file,mcpc_file,startTimesMat);
    
    %Using matchBin.m, combine all xcoords, ycoords, measurements and times
    %into a single matrix (matched_data_mat) and store references to which
    %of those indices correspond to which bin, whose coordinate is
    %represented in central_coord. 
    [W,matched_data_mat,central_coord,count] = matchBin(numBins,time_nmea,xcoord,ycoord,time_mcpc,aveconc,startTimesMat);
    
    %Take the raw data in matched_data_mat and obtain a median value for
    %each bin (accessed using the indices in W). data is a 2D matrix
    %containing a pollution value for each spatial bin (columns) and time (rows).  
    data = ones(size(W,1),size(W,2));
    for i=1:size(W,1)
        for j=1:size(W,2)
            data(i,j) = median(matched_data_mat(W{i,j},4)); 
        end
    end
    
    %log transform data to make it more normal
    data=log(data);
    %transpose data so that time is columns and coordinate is rows.
    transpose_data = data';
    
    %detrend the first two runs (the training sets). Detrend using all data
    %because we assume that all data will be availble to us for the
    %training runs. 
    allTrends = zeros(size(transpose_data,1),size(transpose_data,2));
    transpose_data_detrend = zeros(size(transpose_data,1),size(transpose_data,2));
    for i=1:2
        [allTrends(:,i),transpose_data_detrend(:,i)]=detrend(central_coord',transpose_data(:,i)',true);
    end
    
    %Using constructBasis1D.m, construct Gaussian basis functions. Store
    %these basis functions in S and their coefficients in eta. Coefficients
    %are determined based off of the second training set
    %(transpose_data(:,2)). 
    [S,eta,resolution] = constructBasis1D(central_coord,'Gaussian',r,resolution,transpose_data(:,2));

    %%%%% SET SIGNAL TO NOISE RATIO HERE %%%%%
    SNR=10;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Rename first two data sets as initial_meas_1 and 2. 
    initial_meas_1 = transpose_data_detrend(:,1);
    initial_meas_2 = transpose_data_detrend(:,2);

    %Calculate initial difference (between first two training sets) 
    %normalize by the median (normalized_initial_diff currently not used)
    initial_diff = initial_meas_2-initial_meas_1;
    normalized_initial_diff = initial_diff/median(initial_diff); 
    
    %Remove the training sets from transpose_data. 
    transpose_data=transpose_data(:,3:size(transpose_data,2));
    
    %T, the number of time bins, does NOT include the two parameter estimation runs
    %N is the number of bins
    N=size(transpose_data,1);
    T=size(transpose_data,2); 
    
    %MLvec = missing location vector = the location where you want to start
    %removing data. Initialize as a random location.
    MLvec = 1:T;
    for t=1:T
        MLvec(t)=randi([1 N]);
    end
    %MWvec = missing width vector = the number of bins to be removed.
    %Initialize as zero.
    MWvec = 0.*(1:T);
    
%%%%%%%%%%%     FIT TESTS       %%%%%%%%%%%
    %Used for testing prediction performance 
    
    %Choose the time period to have the data removed from
    missing_one = T; %1 to T
    %Choose where the data removal happens
    missing_loc = 0.3; %0.03 to 1-missing
    %Set all locations in MLvec to be >N (ie, nothing gets removed)
    MLvec(1:T) = N+1;
    %Assign the spatial bin where removal starts in the time bin where you
    %want removal to happen
    MLvec(missing_one) = floor(missing_loc*N);
    %Assign the width of spatiall bins you want removed. 
    MWvec(MLvec<N) = floor(missing*N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Compute the variance-covariance matrix of 1)first training set, 2)
    %the second training set and 3) the lag variance-covariance matrix
    %between the first and second training sets. 
    %NOTE: Because we have already detrended the data (by subtracting off
    %the mean), we can compute variance simply by multiplying by the
    %transpose and dividing by N. 
    sigma_1 = (initial_meas_1*initial_meas_1')/N; 
    sigma_2 = (initial_meas_2*initial_meas_2')/N;
    sigmaLag_2 = ((initial_meas_1*initial_meas_2')/N)';

    %Compute the Q,R decomposition of S (the matrix of spatial basis
    %functions)
    [Q,R]=qr(S,0);
     
    %Calculate Kprelim, the preliminary estimate of K that does not
    %include any information about the errors (D).
    Kprelim = inv(R)*Q'*sigma_1*Q*(inv(R))';
    %Use Kprelim,S and SNR to compute the variances
    sigma2_xi=sum(diag(S*Kprelim*S'))/N;
    sigma2_eps=(1/SNR)*(sum(diag(S*Kprelim*S'))+N*sigma2_xi)/N;
    %By homoscedastic, we mean that variance is not a function of any of
    %our other explanatory variables, so implement this by simply making
    %all variances the same. 
    V_xi_homoscedastic = eye(N); %(this gives you the identity matrix)
    V_xi_heteroscedastic = diag(normalized_initial_diff); %I tried including this but I'm getting Ptt is singular or badly scaled.
    V_eps = eye(N); %(this gives you the identity matrix) 
    %D is the identity matrix multiplied by the sum of both variances
    %(sigma2_xi and sigma2_eps). 
    D= sigma2_xi*V_xi_homoscedastic + sigma2_eps*V_eps;
    
    % sigma_1 is my covariance matrix. D is a matrix representing
    % uncertainty in my spatial process and the measurements from my sensor
    % and eye(N) is the NxN identity matrix. 
    A = D^(-1/2)*sigma_1*D^(-1/2) - eye(N);
    % this line calculates the eigenvectors (ev) and eigenvalues (eigval)
    % of A
    [ev,eigval] = eig(A);
    % this line simply converts my eigenvalue matrix to a 1 vector for
    % simple processing
    d=diag(eigval); 
    % any eigenvalues of A that are less than or equal to 0.1 (arbitrarily
    % chosen by trial and error) are made 0.1
    d(d<=0.1)=0.1; %ISSUE: This is an arbitrary value
    % this line converts new eigenvalues back from vector to matrix
    eigval_PD = diag(d);  
    % reconstruct the new A matrix (A_PD), now positive definite
    A_PD = ev*eigval_PD*ev';
    % reconstruct a new positive definite covariance matrix (sigma_1_PD)
    % from the positive definite A
    sigma_1_PD = (D^(1/2))*A_PD*(D^(1/2)) + D;
    
    %Calculate K_1, the first covariance matrix of eta, the basis function
    %coefficients.
    K_1 = inv(R)*(Q')*(sigma_1_PD-D)*Q*(inv(R)');
    %Calculate L_2, the lag covariance matrix of eta_1 to eta_2, the basis
    %function coefficients for the first and second preliminary runs
    %respectively.
    L_2=inv(R)*(Q')*(sigmaLag_2')*Q*(inv(R)');
    %Calculate the first order autoregressive matrix H, which tells you how
    %eta_t changes to eta_t+1. 
    H_2=L_2*inv(K_1);
    
    %Repeat what you did above but with the covariance matrix for the
    %second run of data (sigma_2). 
    A = D^(-1/2)*sigma_2*D^(-1/2) - eye(N);
    [ev,eigval] = eig(A);
    d=diag(eigval); 
    d(d<=0.1)=0.1; %ISSUE: This is an arbitrary value
    eigval_PD = diag(d); 
    A_PD = ev*eigval_PD*ev';
    sigma_2_PD = (D^(1/2))*A_PD*(D^(1/2)) + D;

    %Calculate K_2, the first covariance matrix of eta, the basis function
    %coefficients for the second preliminary run.
    K_2 = inv(R)*(Q')*(sigma_2_PD-D)*Q*(inv(R)');
    %Calculate U, the covariance matrix for the innovation vector (don't
    %really worry about what that is). 
    U_2=K_2 - H_2*L_2; 

    %We are using static parameters so fix the K, H and U parameters (L
    %gets fixed too, but doesn't get used again). 
    K = K_1; %rxr
    H = H_2; %rxr
    U = U_2; %rxr
    
    %The total number of basis functions is calculated from the resolution
    %and the number of basis functions in the first resolution, r. 
    numBasisFunctions = r*(1-2^resolution)/(-1);
    %reassign r as the TOTAL number of basis functions. Sorry if this is
    %confusing.
    r=numBasisFunctions;
    
    %Initialize empty matrices
    eta_t_t=zeros([r T]); %All basis function coefficient predictions
    eta_tp_t = zeros([r T]); %All basis function coefficient EVOLVED STATES (these are
    %NOT the final predictions, but the intermediate calculation used in
    %the kalman filter). 
    Ptt=zeros([r r T]); %Mean squared errors for eta_t_t
    Ptpt=zeros([r r T]); %Mean squared errors for eta_tp_t
    Y_pred = zeros([N T]); %All predictions of the CONCENTRATIONS
    var_pred = zeros([N T]); %All respective prediction errors for concs
    diff = zeros([N T]); %All differences between measured and pred.
    
    Dii = sigma2_xi+sigma2_eps; %Because D has the same value along the diagonal 
    %we can simplify calcs later on by just assigning that scalar value to
    %Dii (ie, any element along the diagonal).

    %% Kalman Filtering
    %IMPORTANT: This code is not "real time", meaning that it is meant to
    %analyze model performance AFTER you have taken all your data, not to
    %predict based on one time period of data, though it should be very
    %simple to change it to do that. 
    for t=1:T %Loop through all time periods of data you have
        ML = MLvec(t); %ML is the value of MLvec at time t
        MW = MWvec(t); %MW is the value of MWvec at time t
        if MW == 0;
            fullSet = true; %If missing width = 0 (ie, no data is being removed, then fullSet = true.
        else
            fullSet = false;
        end
        
        %Use the detrend function to detrend the data from time t 
        %IMPORTANT: if fullSet is true, then our trend that we subtract off
        %will be based on all the data. If it is false, then only those
        %data that are available (ie, not removed) will be used to
        %calculate the trend (ie, otherwise we would be pretending we had
        %access to information that we really don't). 
        [allTrends(:,t+2),transpose_data_detrend(:,t+2)]=detrend(central_coord',transpose_data(:,t)',fullSet,ML,(ML+MW-1));
        %The sample function will take in the detrended data, ML, and MW and return
        %only those data points that weren't removed (ie, Z_observed). It
        %will also return a sparser representation of S, which has the same
        %number of rows as Z_obs (ie, if the full set of both have N rows, then
        %these have N-MW rows corresponding to the data points that weren't
        %removed). 
        [n, obs_index, Z_obs, St, M] = Sample(transpose_data_detrend(:,t+2),S,ML,MW);
        %for simplicity, reassign Z_obs as Z_t. 
        Z_t = Z_obs;
        %If this is the first time run, then our eta_t_t, eta_tp_t, Ptpt
        %and Ptt matrices are empty so we use eta and K instead of eta_t_t
        %and Ptt respectively. 
        if t==1
            eta_tp_t(:,t)=H*eta; %This is our EVOLVED STATE (ie, how the data should
            %evolve in time from t to t+1 according to our first order
            %autoregressive matrix, H. 
            Ptpt(:,:,t)=H*K*H'+U; %EVOLVED STATE mean squared error
            
            %G=Ptpt(:,:,t)*St'*((1/Dii).*eye(n)-(Dii^(-2)).*St*inv(inv(Ptpt(:,:,t))+(1/Dii).*St'*St)*St');
            %Ptt(:,:,t)=Ptpt(:,:,t)-G*St*Ptpt(:,:,t);
            %eta_t_t(:,t)=eta_tp_t(:,t)+G*(Z_t-St*eta_tp_t(:,t));
            %Calculate 
            %eta_tp_t(:,t)=H*eta_tp_t(:,t);
            %Ptpt(:,:,t)=H*Ptpt(:,:,t)*H'+U;
        else
            %Calculate EVOLVED STATE OF ETA_T using H
            eta_tp_t(:,t)=H*eta_t_t(:,t-1);
            %EVOLVED STATE OF ETA_T mean squared error
            Ptpt(:,:,t)=H*Ptt(:,:,t-1)*H'+U;
        end
            %Calculate Kalman gain matrix, G based on relative errors of
            %MEASUREMENTS (Dii) and EVOLVED STATE OF ETA_T (Ptpt)
            G=Ptpt(:,:,t)*St'*((1/Dii).*eye(n)-(Dii^(-2)).*St*inv(inv(Ptpt(:,:,t))+(1/Dii).*St'*St)*St');
            %Calculate new PREDICTED STATE OF ETA_Tmean squared error using G.
            Ptt(:,:,t)=Ptpt(:,:,t)-G*St*Ptpt(:,:,t);
            %Calculate new PREDICTED STATE OF ETA_T using G
            eta_t_t(:,t)=eta_tp_t(:,t)+G*(Z_t-St*eta_tp_t(:,t));
            %Calculate DETRENDED PREDICTIONS OF CONCENTRATIONS (we have now exited
            %the Kalamn filter)
            Y_pred(:,t)=S*eta_t_t(:,t)+sigma2_xi*M*...
                ((1/Dii).*eye(n)-(Dii^(-2)).*St*inv(inv(Ptpt(:,:,t))+(1/Dii).*St'*St)*St')*...
                (Z_t-St*eta_tp_t(:,t));
            %Add the trend back in to get the PREDICTION OF CONCENTRATIONS
            Y_pred(:,t)=Y_pred(:,t)+allTrends(:,t+2);
            %Calculate the variance of the PREDICTION OF CONCENTRATIONS
            var_pred(:,t)=reshape(diag(S*Ptpt(:,:,t)*S'+sigma2_xi-(sigma2_xi.*M)*...
                ((1/Dii).*eye(n)-(Dii^(-2)).*St*inv(inv(Ptpt(:,:,t))+(1/Dii).*St'*St)*St')*(sigma2_xi.*M')...
                -2*S*Ptpt(:,:,t)*St'*((1/Dii).*eye(n)-(Dii^(-2)).*St*inv(inv(Ptpt(:,:,t))+(1/Dii).*St'*St)*St')*...
                (sigma2_xi.*M')),N,1);     
    end
    %Now we have finished predicting log transformed concentration values
    %for all time periods.
    
    %Remember when we log transformed the data to make it more normal?
    %Well, now we reverse that log transformation (for all time periods we
    %predicted on). 
    Y_pred = exp(Y_pred);
    %Also transform the binned data back.
    transpose_data=exp(transpose_data);
    %Calculate diff (again, this is for all time periods we predicted on).
    diff = Y_pred-transpose_data;
    
    
    %This is strictly for visualizing the predictions. Calculate the range
    %of missing data. 
    missRange = [MLvec(missing_one) (MLvec(missing_one)+MWvec(missing_one)-1)];
    %This is strictly for model testing. Calculate the root mean squared
    %error in the region where data were removed. 
    RMSE_missingReg = sqrt(mean(diff(missRange(1):missRange(2),missing_one).^2));
    
    %This is strictly for visualizing the predictions. Plot all the
    %predictions along with the measured data to see how well our
    %predictions matched the measured data. 
    for i=1:T
        test_case=i;
        availTC=zeros([N 1]);
        availTC(MLvec(test_case):(MLvec(test_case)+MWvec(test_case)))=max(transpose_data(:,test_case));
        figure
        plot(transpose_data(:,test_case));
        %actual is blue
        hold on
        plot(Y_pred(:,test_case));
        %predicted is orange
        plot(availTC(1:N));
        plotTitle = strcat('No Outlier Model Performance on',{' '},startTimesMat(missing_one+2,:));
        title(plotTitle)
        xlabel('Coordinate Bin Number')
        ylabel('Median PM Concentration (counts/cm^{3})')
        legend('Measured Concentrations','Predicted Concentrations');
    end



end