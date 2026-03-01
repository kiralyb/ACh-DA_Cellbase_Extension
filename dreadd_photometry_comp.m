function dreadd_photometry_comp(animalIDlist,ach_mask,da_mask)

sr = 12048;
Win_ACh = sr*2 : sr*2 + sr*0.75;
Win_DA = sr*2 : sr*2 + sr*2;
baselineinx = 1:sr*2;

for i = 1:length(animalIDlist)
    % extract baseline df/f
    load([getpref('cellbase','datapath'),'\',animalIDlist{i},'\',animalIDlist{i},'.mat']);
    baseline_mu_da = mean(mean(F_SO_D{1}(type_SO{1}==0,baselineinx)));
    baseline_std_da = std(mean(F_SO_D{1}(type_SO{1}==0,baselineinx)));
    baseline_mu_ach = mean(mean(F_SO_A{1}(type_SO{1}==0,baselineinx)));
    baseline_std_ach = std(mean(F_SO_A{1}(type_SO{1}==0,baselineinx)));
    for p = 1:3       
        if mod(p,2) == 0
            invfact = -1;
        else
            invfact = 1;
        end
        type_SO_filt = [nan(1,3),type_SO{p}(4:end)];
        type_DF_filt = [nan(1,3),type_DF{p}(4:end)];
        
        for isdread = 1:2
            
            M_SO_A(i,isdread,p,:) = nanmean(((F_SO_A{p}(type_SO_filt==(isdread-1),:))-baseline_mu_ach)./baseline_std_ach);
            M_SO_D(i,isdread,p,:) = nanmean(((F_SO_D{p}(type_SO_filt==(isdread-1),:))-baseline_mu_da)./baseline_std_da);
            [MAX_SO_A(i,isdread,p)] = Photoresponse(M_SO_A(i,isdread,p,Win_ACh),5,0.1,1,sr);
            [MAX_SO_D(i,isdread,p)] = Photoresponse(M_SO_D(i,isdread,p,Win_DA),5,0.1,1,sr);
            
            
            M_DF_A(i,isdread,p,:) = nanmean(((F_DF_A{p}(type_DF_filt==(isdread-1),:))-baseline_mu_ach)./baseline_std_ach);
            M_DF_D(i,isdread,p,:) = nanmean(((F_DF_D{p}(type_DF_filt==(isdread-1),:))-baseline_mu_da)./baseline_std_da);
            [MAX_DF_A(i,isdread,p)] = Photoresponse(M_DF_A(i,isdread,p,Win_ACh),5,0.1,1,sr);
            [MAX_DF_D(i,isdread,p)] = Photoresponse(M_DF_D(i,isdread,p,Win_DA),5,0.1,invfact,sr);
            
            M_DF_A_abs(i,isdread,p,:) = nanmean(((F_DF_A{p}(type_DF_filt==(isdread-1),:))));%-baseline_mu_ach)./baseline_std_ach);
            M_DF_D_abs(i,isdread,p,:) = nanmean(((F_DF_D{p}(type_DF_filt==(isdread-1),:))));%-baseline_mu_da)./baseline_std_da);
            [~,MAX_DF_A_abs(i,isdread,p)] = Photoresponse(M_DF_A_abs(i,isdread,p,Win_ACh),5,0.1,1,sr);
            [~,MAX_DF_D_abs(i,isdread,p)] = Photoresponse(M_DF_D_abs(i,isdread,p,Win_DA),5,0.1,invfact,sr);
            
        end
        
    end
    
end

% 4f
DREADD_Comp(M_SO_A,MAX_SO_A,ach_mask,3)
xlabel(['Time from ','Cue',' (s)'])
% 4h
DREADD_Comp(M_SO_D,MAX_SO_D,da_mask,3)
xlabel(['Time from ','Cue',' (s)'])
% 4i
C21_COMP(squeeze(M_DF_D_abs(:,2,:,:)),squeeze(MAX_DF_D_abs(:,2,:)),da_mask,2)
xlabel(['Time from ','Feedback',' (s)'])
% 4j
DREADD_Comp(M_DF_D,MAX_DF_D,da_mask,2)
xlabel(['Time from ','Feedback',' (s)'])

% S10g
C21_COMP(squeeze(M_DF_A_abs(:,2,:,:)),squeeze(MAX_DF_A_abs(:,2,:)),ach_mask,2)
xlabel(['Time from ','Feedback',' (s)'])
DREADD_Comp(M_DF_A,MAX_DF_A,ach_mask,2)
xlabel(['Time from ','Feedback',' (s)'])
% S10i
DREADD_Comp(M_DF_A,MAX_DF_A,ach_mask,3)
xlabel(['Time from ','Feedback',' (s)'])
DREADD_Comp(M_DF_D,MAX_DF_D,da_mask,3)
xlabel(['Time from ','Feedback',' (s)'])
