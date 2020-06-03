clear; close all

%% ECHAM LIG

load('LIG.mat')
load('LIG_PI.mat')

wiso_anom = wiso_LIG - wiso_LIG_PI;
tot_pre_anom = pre_tot_ann_LIG - pre_tot_ann_LIG_PI;
pre_anom = pre_ann_LIG - pre_ann_LIG_PI;
temp_anom = temp_LIG - temp_LIG_PI;

a = wiso_anom(49:96,:); b = wiso_anom(1:48,:); wiso_anom = [a;b];
clear a b; a = pre_anom(49:96,:); b = pre_anom(1:48,:); pre_anom = [a;b];
clear a b; a = tot_pre_anom(49:96,:); b = tot_pre_anom(1:48,:); tot_pre_anom = [a;b];
clear a b; a = temp_anom(49:96,:); b = temp_anom(1:48,:); temp_anom = [a;b];

matrix_resized = imresize(wiso_anom, [320 160], 'nearest'); LIG_wiso_anom = matrix_resized; clear matr*
matrix_resized = imresize(pre_anom, [320 160], 'nearest'); LIG_pre_anom = matrix_resized; clear matr*
matrix_resized = imresize(tot_pre_anom, [320 160], 'nearest'); LIG_tot_pre_anom = matrix_resized; clear matr*
matrix_resized = imresize(temp_anom, [320 160], 'nearest'); LIG_temp_anom = matrix_resized; clear matr*

load('mask_LIG.mat')
a = slm(49:96,:); b = slm(1:48,:); slm = [a;b];
matrix_resized = imresize(slm, [320 160], 'nearest'); land_mask_LIG = matrix_resized; clear matr*

save('ECHAM_dat','LIG_wiso_anom','LIG_pre_anom','LIG_tot_pre_anom','LIG_temp_anom','land_mask_LIG')


%% ECHAM MH and LGM
clear

load('MH_files.mat'); load('MH_PI_files.mat')
MH_wiso_anom = wiso_MH_ave - wiso_MH_PI_ave;
MH_pre_anom = pre_MH_ave - pre_MH_PI_ave;
MH_tot_pre_anom = pre_tot_MH_ave - pre_tot_MH_PI_ave;
MH_temp_anom = temp_MH_ave - temp_MH_PI_ave;

a = MH_wiso_anom(161:320,:); b = MH_wiso_anom(1:160,:); MH_wiso_anom = [a;b];
clear a b; a = MH_pre_anom(161:320,:); b = MH_pre_anom(1:160,:); MH_pre_anom = [a;b];
clear a b; a = MH_tot_pre_anom(161:320,:); b = MH_tot_pre_anom(1:160,:); MH_tot_pre_anom = [a;b];
clear a b; a = MH_temp_anom(161:320,:); b = MH_temp_anom(1:160,:); MH_temp_anom = [a;b];

save('ECHAM_dat','MH_wiso_anom','MH_pre_anom','MH_tot_pre_anom','MH_temp_anom','-append')


load('LGM_files.mat'); load('LGM_PI_files.mat')
LGM_wiso_anom = wiso_LGM_ave - wiso_LGM_PI_ave;
LGM_pre_anom = pre_LGM_ave - pre_LGM_PI_ave;
LGM_tot_pre_anom = pre_tot_LGM_ave - pre_tot_LGM_PI_ave;
LGM_temp_anom = temp_LGM_ave - temp_LGM_PI_ave;

a = LGM_wiso_anom(161:320,:); b = LGM_wiso_anom(1:160,:); LGM_wiso_anom = [a;b];
clear a b; a = LGM_pre_anom(161:320,:); b = LGM_pre_anom(1:160,:); LGM_pre_anom = [a;b];
clear a b; a = LGM_tot_pre_anom(161:320,:); b = LGM_tot_pre_anom(1:160,:); LGM_tot_pre_anom = [a;b];
clear a b; a = LGM_temp_anom(161:320,:); b = LGM_temp_anom(1:160,:); LGM_temp_anom = [a;b];


save('ECHAM_dat','LGM_wiso_anom','LGM_pre_anom','LGM_tot_pre_anom','LGM_temp_anom','-append')
