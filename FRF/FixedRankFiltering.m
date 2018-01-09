function [avr,avm,avb,Y_pred,var_pred,diff,transpose_data] = FixedRankFiltering(nmea_file,mcpc_file,startTimesMat,missing,r,numBins)
%Parameters to test: The time of run, the distribution of the data (ie, are the data points skewed,normal? and does this differ from run to run?)
%More parameters to test: number of bins, number of missing data points,
%location of missing data points (ie, was it obscuring a spike), rank, number of resolutions, SNR

%possible measurements of goodness: diff, diff at locations of missing
%data,var_pred,avr,avm,avb,(and all three of these at locations of missing 
%data), direction of difference OTHERS

%Should also be testing how well the binning works, and how well the
%detrending works (could use leave on out method or cross-validation could 
%work for testing how well the linear regression worked in the detrend
%function).

%Make sure to give credit to Kuztfuss and Cressie for using some of their
%code, as well as Cressie et al for using FRF algorithm and part of their code
%Relies on helper functions Sample.m, constructBasis1D.m, readData.m,
%matchBin.m and DetermineRuns.m
%Currently only works for 1D. And assumes the preliminary runs contain all data for all the bins.
%missing is the percentage of observations that are excluded in each time period from the full set of observations
%r is the number of basis functions to use relative to the number of grid  locations (expressed as a percentage)
    
    [xcoord,ycoord,time_nmea,aveconc,time_mcpc] = readData(nmea_file,mcpc_file,startTimesMat);
    [W,data_detrend,X0,betahat,central_coord,count] = matchBin(numBins,time_nmea,xcoord,ycoord,time_mcpc,aveconc,startTimesMat);
    data = ones(size(W,1),size(W,2));
    trend = X0*betahat;
    large_scale_trend = data;
    for i=1:size(W,1)
        for j=1:size(W,2)
            data(i,j) = mean(data_detrend(W{i,j},5));
            large_scale_trend(i,j) = mean(trend(W{i,j}));
        end
    end
    %data contains the DETRENDED aveconc. Detrending took place with full
    %knowledge of all the data
    %I'll need to add the trend back in at the end. 
    
    %transpose data so that time is rows and coordinate is columns.
    transpose_data = data';
    
    [S,eta] = constructBasis1D(central_coord,'Gaussian',r,1,transpose_data(:,1));
    
    SNR=5;
    initial_meas_1 = transpose_data(:,1);
    initial_meas_2 = transpose_data(:,2);
    N=size(transpose_data,1);
    T=size(transpose_data,2);
    MLvec = 1:T;
    MLvec(MLvec>0)=N+1;
    for t=1:T
        MLvec(t)=randi([1 N]);
    end
    MWvec = 0.*(1:T);
    MWvec(MLvec<N) = floor(missing*N); 
    sigma_1 = (initial_meas_1'*initial_meas_1)/N; %This computes covariance because data is already detrended 
    sigma_2 = (initial_meas_2'*initial_meas_2)/N;
    sigmaLag_2 = ((initial_meas_1'*initial_meas_2)/N)';
    
    data=data(3:T,:);

    [Q,R]=qr(S,0);
    Kprelim = inv(R)*Q'*sigma_1*Q*(inv(R))';
    %Should make sigma2_eps have to do with SNR
    %what should sigma_xi be?
    sigma2_xi=sum(diag(S*Kprelim*S'))/N;
    sigma2_eps=(1/SNR)*(sum(diag(S*Kprelim*S'))+N*sigma2_xi)/N;
    V_xi = eye(N); %assumes no heteroskedacticity
    V_eps = eye(N);
    D= sigma2_xi*V_xi + sigma2_eps*V_eps;
    
    A = D^(-1/2)*sigma_1*D^(-1/2) - eye(N);
    [ev,eigval] = eig(A);
    d=diag(eigval); 
    d(d<=0.1)=0.1; %ISSUE: This is an arbitrary value
    eigval_PD = diag(d); 
    A_PD = ev*eigval_PD*ev';
    sigma_1_PD = (D^(1/2))*A_PD*(D^(1/2)) + D;
    
    K_1 = inv(R)*(Q')*(sigma_1_PD-D)*Q*(inv(R)');
    L_2=inv(R)*(Q')*(sigmaLag_2')*Q*(inv(R)');
    H_2=L_2*inv(K_1);
    
    A = D^(-1/2)*sigma_2*D^(-1/2) - eye(N);
    [ev,eigval] = eig(A);
    d=diag(eigval); 
    d(d<=0.1)=0.1; 
    eigval_PD = diag(d); 
    A_PD = ev*eigval_PD*ev';
    sigma_2_PD = (D^(1/2))*A_PD*(D^(1/2)) + D;

    K_2 = inv(R)*(Q')*(sigma_2_PD-D)*Q*(inv(R)');
    U_2=K_2 - H_2*L_2; 

    K = K_1; %rxr
    H = H_2; %rxr
    U = U_2; %rxr

    eta_t_t=zeros([r T]);
    eta_t_t(:,1)=eta;
    Y_pred = zeros([N T]);
    Y_pred(:,1) = S*eta_t_t(:,1);
    var_pred = zeros([N T]);
    diff = zeros([N T]);
    diff(:,1) = Y_pred(:,1)-transpose_data(:,1);
    Dii = sigma2_xi+sigma2_eps;

    Ptt=zeros([r r T]);
    Ptt(:,:,1)=K;
    Ptpt=zeros([r r T]);
    eta_tp_t = zeros([r T]);
    
    for t=1:T
        ML = MLvec(t);
        MW = MWvec(t);
        [n, obs_index, Z_obs, St, M] = Sample(transpose_data(:,t),S,ML,MW);
        Z_t = Z_obs;
        
        if t==1
            eta_tp_t(:,t)=H*eta_t_t(:,1);
            Ptpt(:,:,t)=H*K*H'+U;
            for i=1:3 %Allow for convergence
                G=Ptpt(:,:,t)*St'*((1/Dii).*eye(n)-(Dii^(-2)).*St*inv(inv(Ptpt(:,:,t))+(1/Dii).*St'*St)*St');
                Ptt(:,:,t)=Ptpt(:,:,t)-G*St*Ptpt(:,:,t);
                eta_t_t(:,t)=eta_tp_t(:,t)+G*(Z_t-St*eta_tp_t(:,t));
                eta_tp_t(:,t)=H*eta_tp_t(:,t);
                Ptpt(:,:,t)=H*Ptpt(:,:,t)*H'+U;
            end
        else
            eta_tp_t(:,t)=H*eta_t_t(:,t-1);
            Ptpt(:,:,t)=H*Ptt(:,:,t-1)*H'+U;
        end
            G=Ptpt(:,:,t)*St'*((1/Dii).*eye(n)-(Dii^(-2)).*St*inv(inv(Ptpt(:,:,t))+(1/Dii).*St'*St)*St');
            Ptt(:,:,t)=Ptpt(:,:,t)-G*St*Ptpt(:,:,t);
            eta_t_t(:,t)=eta_tp_t(:,t)+G*(Z_t-St*eta_tp_t(:,t));
            Y_pred(:,t)=S*eta_t_t(:,t)+sigma2_xi*M*...
                ((1/Dii).*eye(n)-(Dii^(-2)).*St*inv(inv(Ptpt(:,:,t))+(1/Dii).*St'*St)*St')*...
                (Z_t-St*eta_tp_t(:,t));
            var_pred(:,t)=reshape(diag(S*Ptpt(:,:,t)*S'+sigma2_xi-(sigma2_xi.*M)*...
                ((1/Dii).*eye(n)-(Dii^(-2)).*St*inv(inv(Ptpt(:,:,t))+(1/Dii).*St'*St)*St')*(sigma2_xi.*M')...
                -2*S*Ptpt(:,:,t)*St'*((1/Dii).*eye(n)-(Dii^(-2)).*St*inv(inv(Ptpt(:,:,t))+(1/Dii).*St'*St)*St')*...
                (sigma2_xi.*M')),N,1);     
            diff(:,t)=Y_pred(:,t)-transpose_data(:,t);
    end
    Y_pred = (Y_pred'+large_scale_trend)';
    transpose_data = (transpose_data'+large_scale_trend)';
    
    [r,m,b]=regression(Y_pred',transpose_data');
    avr = mean(r);
    avm = mean(m);
    avb = mean(b);
    
    if T>10
        warning('You are trying to plot more than 10 things')
    end
    
    for i=1:T
        test_case=i;
        availTC=zeros([N 1]);
        availTC(MLvec(test_case):(MLvec(test_case)+MWvec(test_case)))=max(data(:,test_case));
        figure
        plot(Y_pred(:,test_case));
        hold on
        plot(transpose_data(:,test_case));
        plot(availTC(1:N));
    end
    
    figure
    surf(diff)
    
    figure
    surf(transpose_data)
    
end