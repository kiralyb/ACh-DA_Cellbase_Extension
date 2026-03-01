function firstnew_photometry_run(animalIDlist,ach_mask,da_mask)

sr = 12048;
Win_ACh = sr*2 : sr*2 + sr*0.75;
Win_DA = sr*2 : sr*2 + sr*2;

for i = 1:length(animalIDlist) 
    load([getpref('cellbase','datapath'),'\',animalIDlist{i},'\',animalIDlist{i},'.mat']);
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