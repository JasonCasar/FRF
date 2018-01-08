\\\function [avr,avm,avb,Y_pred,var_pred,diff] = FixedRankFiltering(nmea_file,mcpc_file,startTimesMat,missing,r,numBins)
%Make sure to give credit to Kuztfuss and Cressie for using some of their
%code, as well as Cressie et al for using FRF algorithm and part of their code
%Relies on helper functions Sample.m, constructBasis1D.m, readData.m,
%matchBin.m and DetermineRuns.m
%Currently only works for 1D. And assumes the preliminary runs contain all data for all the bins.
%missing is the percentage of observations that are excluded in each time period from the full set of observations
%r is the number of basis functions to use relative to the number of grid  locations (expressed as a percentage)
    
    [xcoord,ycoord,time_nmea,aveconc,time_mcpc] = readData(nmea_file,mcpc_file,startTimesMat);
    [W,data_detrend,X0,betahat,central_coord,count] = matchBin(numBins,time_nmea,xcoord,ycoord,time_mcpc,aveconc,startTimesMat);
    %data_detrend(W{1,1},5)
    
    %In the matlab code for the fixed rank filtering tutorial, they detrend
    %their data before binning it (they did their regression on longitude. 
    %see lines 25-32 of Xco2analysis), which means their covariance matrix 
    %is computed by averaging each of the residuals in each bin, then
    %squaring to get the variance (if the same bin) and multiplying both
    %to get the covariance (if separate bins). I should write a function to
    %detrend the data, and return a cell of matrices. I will need to look at the xco2analysis functions to get a good idea of how to do this.  
    
    SNR=10;
    initial_meas_1 = data(1,:);
    initial_meas_2 = data(1,:);
    N=length(initial_meas_1);
    T=length(data(:,1));
    MLvec = 1:T;
    MLvec(MLvec>0)=N+1;
    for t=1:T
        MLvec(t)=randi([1 N]);
    end
    MWvec = 0.*(1:T);
    MWvec(MLvec<N) = floor(missing*N); 
    sigma_1 = cov(initial_meas_1,1)';
    sigma_2 = cov(initial_meas_2,1)';
    sigmaLag_2 = ((1/N)*(initial_meas_1-(1/N)*eye(obs)*initial_meas_1)'*(initial_meas_2-(1/N)*eye(obs)*initial_meas_2))';
    
    data=data(3:T,:);
    [S,eta] = constructBasis1D(1:N,'Gaussian',r,1,(mean(initial_meas_2))');

    [Q,R]=qr(S,0);
    Kprelim = inv(R)*Q'*sigma_1*Q*(inv(R))';
    %Should make sigma2_eps have to do with SNR
    %Should make sigma2_xi be variance of pollution value in each cell?
    sigma2_xi=sum(diag(S*Kprelim*S'))/N;
    sigma2_eps=(1/SNR)*(sum(diag(S*Kprelim*S'))+N*sigma2_xi)/N;
    D= sigma2_xi*eye(N) + sigma2_eps*eye(N); %assumes no heteroskedacticity
    
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
    diff(:,1) = Y_pred(:,1)-data(:,1);
    Dii = sigma2_xi+sigma2_eps;

    Ptt=zeros([r r T]);
    Ptt(:,:,1)=K;
    Ptpt=zeros([r r T]);
    eta_tp_t = zeros([r T]);
    
    for t=1:T
        ML = MLvec(t);
        MW = MWvec(t);
        [n, obs_index, Z_obs, St, M] = Sample(data(:,t),S,ML,MW);
        Z_t = Z_obs;
        
        if t==1
            eta_tp_t(:,t)=H*coeffs;
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
            diff(:,t)=Y_pred(:,t)-data(:,t);
    end
    
    [r,m,b]=regression(Y_pred',data');
    avr = mean(r);
    avm = mean(m);
    avb = mean(b);
    
    if N>10
        warning('You are trying to plot more than 10 things')
    end
    
    for i=1:T
        test_case=i;
        availTC=zeros([N 1]);
        availTC(MLvec(test_case):(MLvec(test_case)+MWvec(test_case)))=max(data(:,test_case));
        figure
        plot(Y_pred(:,test_case));
        hold on
        plot(data(:,test_case));
        plot(availTC(1:N));
    end
    
    figure
    surf(diff)
    
end