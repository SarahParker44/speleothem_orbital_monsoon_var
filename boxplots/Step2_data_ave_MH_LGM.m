%This code produces the MH_files.mat, LGM_files.mat files with averages during those
%periods for each of the variables (total annual precipitation, nanmean
%surface temperature, precipitation and precipitation weighted d18O

clear 
close all

for z=1:4

        if z==1
            clear
            load('MH.mat')
            pre_MH_ave=nanmean(pre_ann_MH,3);
            pre_tot_MH_ave=nanmean(pre_tot_ann_MH,3);
            temp_MH_ave=nanmean(temp_MH,3);
            wiso_MH_ave=nanmean(wiso_MH,3);
            save('MH_files.mat')
            
        elseif z==2
            clear
            load('MH_PI.mat')
            pre_MH_PI_ave=nanmean(pre_ann_MH_PI,3);
            pre_tot_MH_PI_ave=nanmean(pre_tot_ann_MH_PI,3);
            temp_MH_PI_ave=nanmean(temp_MH_PI,3);
            wiso_MH_PI_ave=nanmean(wiso_MH_PI,3);
            save ('MH_PI_files.mat')
            
        elseif z==3
            clear
            load('LGM.mat')
            pre_LGM_ave=nanmean(pre_ann_LGM,3);
            pre_tot_LGM_ave=nanmean(pre_tot_ann_LGM,3);
            temp_LGM_ave=nanmean(temp_LGM,3);
            wiso_LGM_ave=nanmean(wiso_LGM,3);
            save ('LGM_files.mat')
            
        elseif z==4
            clear
            load('LGM_PI.mat')
            pre_LGM_PI_ave=nanmean(pre_ann_LGM_PI,3);
            pre_tot_LGM_PI_ave=nanmean(pre_tot_ann_LGM_PI,3);
            temp_LGM_PI_ave=nanmean(temp_LGM_PI,3);
            wiso_LGM_PI_ave=nanmean(wiso_LGM_PI,3);
            save ('LGM_PI_files.mat')
        end  
 
end

disp('You can now delete all files in this folder that are not named as "XX_files.mat"')