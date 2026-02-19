%% loading
choosecb('Cellbase_psychometric')
cellbasepath = getpref('cellbase','datapath')
loadcb
animalIDlist = CELLIDLIST;
type_animal = TheMatrix;
n=length(animalIDlist);
sr = 12048;
baselineinx = 1:sr*2;
ach_mask = [0,0,1,1,1,nan,0,0,0,1,1,1,0,nan];
Win_ACh = sr*2 : sr*2 + sr*0.75;
da_mask = [0,0,1,1,1,nan,0,0,0,1,1,1,nan,0];
Win_DA = sr*2 : sr*2 + sr*2;
%%
novelbehav_DREADD
%%
% 4c
[DIFF_D] = psychometric_DREADD_COMP(animalIDlist(type_animal==1),ones(1,sum(type_animal==1))*4);
% S10c
[DIFF_C] = psychometric_DREADD_COMP(animalIDlist(type_animal==0),ones(1,sum(type_animal==1))*4);
% S10d
boxstat(DIFF_C,DIFF_D,'control','dreadd')
%% S10e
for i = 1:n
    load([cellbasepath,'\',animalIDlist{i},'\',animalIDlist{i},'.mat']);
    % extract C21 (1) and Control (0) days
    MHit(i,1) = mean(Hit(type_SO{1}==0));
    MHit(i,2) = mean(Hit(type_SO{1}==1));
    MFA(i,1) = mean(FalseAlarm(type_SO{1}==0));
    MFA(i,2) = mean(FalseAlarm(type_SO{1}==1));
    MRT(i,1) = mean(RT(type_SO{1}==0));
    MRT(i,2) = mean(RT(type_SO{1}==1));
    
end

fix_DREADD_comp(MRT,type_animal)
ylabel('Reaction time (s)')
fix_DREADD_comp(MHit,type_animal)
ylabel('Hit %')
fix_DREADD_comp(MFA,type_animal)
ylabel('FA %')
%% 
for i = 1:n 
    load([cellbasepath,'\',animalIDlist{i},'\',animalIDlist{i},'.mat']);
    for p = 1:3         
            firsti = find(type_SO{p}==1,1,'first');
            M_SO_A7(i,p,:) = F_SO_A{p}(firsti,:);
            M_SO_D7(i,p,:) = F_SO_D{p}(firsti,:);
            [~,MAX_SO_A7(i,p)] = Photoresponse(M_SO_A7(i,p,Win_ACh),5,0.1,1,sr);
            [~,MAX_SO_D7(i,p)] = Photoresponse(M_SO_D7(i,p,Win_DA),5,0.1,1,sr);          
    end 
end
% Fig 4e, S10f
firstnew_photometry(M_SO_A7,MAX_SO_A7,ach_mask)
% Fig 4g, S10f
firstnew_photometry(M_SO_D7,MAX_SO_D7,da_mask)
%% load
for i = 1:n
    % extract baseline df/f
    load([cellbasepath,'\',animalIDlist{i},'\',animalIDlist{i},'.mat']);
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
%%
% 4f
DREADD_Comp(M_SO_A,MAX_SO_A,ach_mask,3)
% 4h
DREADD_Comp(M_SO_D,MAX_SO_D,da_mask,3)
% 4i
C21_COMP(squeeze(M_DF_D_abs(:,2,:,:)),squeeze(MAX_DF_D_abs(:,2,:)),da_mask,2)
% 4j
DREADD_Comp(M_DF_D,MAX_DF_D,da_mask,2)
% S10g
C21_COMP(squeeze(M_DF_A_abs(:,2,:,:)),squeeze(MAX_DF_A_abs(:,2,:)),ach_mask,2)
DREADD_Comp(M_DF_A,MAX_DF_A,ach_mask,2)
% S10i
DREADD_Comp(M_DF_A,MAX_DF_A,ach_mask,3)
DREADD_Comp(M_DF_D,MAX_DF_D,da_mask,3)

