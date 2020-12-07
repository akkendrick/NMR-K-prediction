function [b] = bootstrap_nmr_called_2 (data, n, m)
% this function is called by bootstrp.m. x is a m x (p+1) matrix, where p is
% the number of predictor variables. The last column of data should be permeability. 
% For NMR data, x should be [log10(T2ML)] or [log10(T2ML), log10(phi)]. 
% If n or m should be fixed to a constant value, then supply them. Otherwise 
% the code will bootstrap over whatever parameters are needed for the 
% supplied x matrix. 

    if size(data, 2) > size(data,1)
        data = data'; 
    end
    if size(data,2) == 3 
        x = [ones(size(data,1),1), data(:,1), data(:,2)]; 
    elseif size(data,2) == 2
        x = [ones(size(data,1),1), data(:,1)]; 
    end

    y = data(:,end);
    
    if nargin == 1
        b = regress(y, x); 
    else
        if nargin == 3
            ft = fittype('poly11');
            cterm1 = m*x(:,3);
            cterm2 = n*x(:,2);
            x = [cterm1, cterm2]; 
            fo = fitoptions('method','LinearLeastSquares','Lower',[-Inf, 1, 1],'Upper',[Inf, 1, 1]);
            model = fit(x,y,ft,fo);
            b(1) = model.p00;  
            b(2) = model.p01;
            b(3) = model.p10;
        elseif nargin == 2
            ft = fittype('poly11');
            cterm1 = x(:,3); 
            cterm2 = n*x(:,2);
            x = [cterm1, cterm2]; 
            %fo = fitoptions('method','LinearLeastSquares','Lower',[-Inf, 1, 1],'Upper',[Inf, 1, 1]);
            %fo = fitoptions('method','LinearLeastSquares','Lower',[-Inf, -2, 1],'Upper',[Inf, 4, 1]);
            fo = fitoptions('method','LinearLeastSquares','Lower',[-Inf, -100, 1],'Upper',[Inf, 100, 1]);

            model = fit(x,y,ft,fo); 
            b(1) = model.p00;
            b(2) = n;
            b(3) = model.p10; 
        else        
            error('Bad x')
        end       
    end
end
