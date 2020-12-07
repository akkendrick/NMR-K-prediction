% Compare parameter distributions


siteNames = {'wisc_all','Site1-WellG5','Site1-WellG6','Site2-WellPN1',...
    'Site2-WellPN2'};


for kk = 1:length(siteNames)
   siteNames{kk}
   computeProfile(siteNames{kk},[],[],0,1);   
end
%%
close all

saveNames = {'wisc_all','G5','G6','PN1','PN2'};
legendNames = {'All','G5','G6','PN1','PN2'};


for kk = 1:length(saveNames)
   name1 = strcat(saveNames{kk},'_bootstrap_n_m_var.mat');
   name2 = strcat(saveNames{kk},'_basic_solving_n_m_var.mat');
   name3 = strcat(saveNames{kk},'_MCMC_n_m_var.mat');
   
   load(name1)
   load(name2)
   load(name3)
   
   figure(1)
   subplot(2,2,1)
   hold on
   histogram(blog_mcmc,'Normalization','pdf')
   xlabel('log_{10}(b)')
   box on
   grid on
      set(gca,'FontSize',14)

   subplot(2,2,2)
   hold on
   histogram(n_mcmc,'Normalization','pdf')
   xlabel('n')
   box on
   grid on
      set(gca,'FontSize',14)

   subplot(2,2,3)
   hold on
   histogram(m_mcmc,'Normalization','pdf')
   xlabel('m')
   box on
   grid on
   legend(legendNames)
   
   set(gca,'FontSize',14)
   
   %%%%%%%%%%%%%%%%%
   logbs = log10(bs_basic);
   
   figure(2)
   subplot(2,2,1)
   hold on
   histogram(logbs,'Normalization','pdf')
   xlabel('log_{10}(b)')
   box on
   grid on
      set(gca,'FontSize',14)

   subplot(2,2,2)
   hold on
   stem(nDirect)
   xlabel('n')
   box on
   grid on
      set(gca,'FontSize',14)

   subplot(2,2,3)
   hold on
   stem(mDirect)
   xlabel('m')
   box on
   grid on
   legend(legendNames)

   set(gca,'FontSize',14)

   
   %%%%%%%%%%%%%%%%%
   figure(3)
   subplot(2,2,1)
   hold on
   histogram(bs,'Normalization','pdf')
   xlabel('log_{10}(b)')
   box on
   grid on
      set(gca,'FontSize',14)

   subplot(2,2,2)
   hold on
   histogram(n_boot,'Normalization','pdf')
   xlabel('n')
   box on
   grid on
      set(gca,'FontSize',14)

   subplot(2,2,3)
   hold on
   histogram(m_boot,'Normalization','pdf')
   xlabel('m')
   box on
   grid on
   legend(legendNames)
   set(gca,'FontSize',14)

   
   
    
    
end