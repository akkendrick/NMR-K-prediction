function  [k_estimates_n1, k_estimates_n2, avgK_n1, avgK_n2] = computeAvgKModel(siteListInd,b,m,n,phi,T2ML)

    SDR_K = @(b,m,n,phi,T2ML) b.*(phi.^m).*(T2ML).^n;

   k_estimates_n1 = [];
   k_estimates_n2 = [];
    Ktemp = {};
    
    for h = 1:length(n)       
        for j = 1:length(siteListInd)
            Ktemp{j} = SDR_K(b(h),m(h),n(h),phi{j},T2ML{j});
        end
        
        if n(h) == 1
            k_estimates_n1 = [k_estimates_n1 vertcat(Ktemp{:})];
        elseif n(h) == 2
            k_estimates_n2 = [k_estimates_n2 vertcat(Ktemp{:})];
        end
        
    end
    

   if ~isempty(k_estimates_n1)
       avgK_n1 = mean(k_estimates_n1,2);
   else
       avgK_n1 = [];
   end
   
   if ~isempty(k_estimates_n2)
      avgK_n2 = mean(k_estimates_n2,2);
   else
       avgK_n2 = [];
   end

    
end

