function [data_detrend,X0,betahat] = detrend(matched_data_mat,startIndices)
%detrend all data given complete knowledge of all data. This also
%eliminates data outside of the start of the first run and the end of the
%last run.

matched_data_mat = matched_data_mat(startIndices(1):startIndices(length(startIndices)),:);
X0 = [ones(size(matched_data_mat,1),1) matched_data_mat(:,2)];
betahat = regress(matched_data_mat(:,4),X0);
data_detrend = horzcat(matched_data_mat,matched_data_mat(:,4)-X0*betahat);

end
