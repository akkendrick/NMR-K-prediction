function totalError = computeError(DPP_K,estimated_K)
% Compute RMSE between model and DPP measurements

   RMSE = sqrt((estimated_K - DPP_K).^2);
   totalError = sum(RMSE);

%    MAE = sum(abs(estimated_K-DPP_K))/length(DPP_K);
%    totalError = MAE;
   
%    MAE = median(abs(estimated_K-DPP_K))/length(DPP_K);
%    totalError = MAE;
end

