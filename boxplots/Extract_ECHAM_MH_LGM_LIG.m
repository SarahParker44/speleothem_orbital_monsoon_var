clear; close all

disp('Extracting d18O for MH runs')

clear path; path=('C:\Users\ph805612\Documents\SISAL\R_programming\Data\Echam5-wiso.T106_EXP011_PL_6k.global.monmean.d18O_prec.nc');
    %ncdisp(path) % check variable name is correct
    isoMH=ncread(path,'wisoaprt_d'); isoMH=squeeze(isoMH);
    time=ncread(path,'time');
    lon=ncread(path,'lon'); lat=ncread(path,'lat');
    
Tbase = datenum(1954, 1, 1);
duration=sscanf(datestr(time(end,1) + Tbase),'15-Dec-%d')-sscanf(datestr(time(1,1) + Tbase),'15-Jan-%d')+1

%% 

clear Tbase duration

clear path; path=('C:\Users\ph805612\Documents\SISAL\R_programming\Data\Echam5-wiso.T106_EXP010_PL_PI.global.monmean.d18O_prec.nc');
    %ncdisp(path) % check variable name is correct
    isoMH_PI=ncread(path,'wisoaprt_d'); isoMH_PI=squeeze(isoMH_PI);

clear path; path=('C:\Users\ph805612\Documents\SISAL\R_programming\Data\Echam5-wiso.T106_EXP011_PL_6k.global.monmean.prec.nc');
    %ncdisp(path) % check variable name is correct
    preMH=ncread(path,'aprt'); preMH=preMH*(2.628*(10^6)); 
    % pre = units in kg/m^2s => convert to mm/month

clear path; path=('C:\Users\ph805612\Documents\SISAL\R_programming\Data\Echam5-wiso.T106_EXP010_PL_PI.global.monmean.prec.nc');
    %ncdisp(path) % check variable name is correct
    preMH_PI=ncread(path,'aprt'); preMH_PI=preMH_PI*(2.628*(10^6)); 
% pre = units in kg/m^2s => convert to mm/month

clear path; path=('C:\Users\ph805612\Documents\SISAL\R_programming\Data\Echam5-wiso.T106_EXP011_PL_6k.global.monmean.tsurf.nc');
    %ncdisp(path) % check variable name is correct
    tmpMH=ncread(path,'tsurf'); 
    
clear path; path=('C:\Users\ph805612\Documents\SISAL\R_programming\Data\Echam5-wiso.T106_EXP010_PL_PI.global.monmean.tsurf.nc');
    %ncdisp(path) % check variable name is correct
    tmpMH_PI=ncread(path,'tsurf');
    
  
    clear path
    save Runs_MH
    
    %% 
    clear; close all

disp('Extracting d18O for LGM runs')

clear path; path=('C:\Users\ph805612\Documents\SISAL\R_programming\Data\Echam5-wiso.T106_EXP016_LGM.global.monmean.d18O_prec.nc');
    %ncdisp(path) % check variable name is correct
    isoLGM=ncread(path,'wisoaprt_d'); isoLGM=squeeze(isoLGM);
    time=ncread(path,'time');
    lon=ncread(path,'lon'); lat=ncread(path,'lat');
 
Tbase = datenum(1978, 1, 1);
duration=sscanf(datestr(time(end,1) + Tbase),'14-Dec-%d')-sscanf(datestr(time(1,1) + Tbase),'15-Jan-%d')+1
clear Tbase duration    

clear path; path=('C:\Users\ph805612\Documents\SISAL\R_programming\Data\Echam5-wiso.T106_EXP015_PD.global.monmean.d18O_prec.nc');
    %ncdisp(path) % check variable name is correct
    isoLGM_PI=ncread(path,'wisoaprt_d'); isoLGM_PI=squeeze(isoLGM_PI);

clear path; path=('C:\Users\ph805612\Documents\SISAL\R_programming\Data\Echam5-wiso.T106_EXP016_LGM.global.monmean.prec.nc');
    %ncdisp(path) % check variable name is correct
    preLGM=ncread(path,'aprt'); 
    
clear path; path=('C:\Users\ph805612\Documents\SISAL\R_programming\Data\Echam5-wiso.T106_EXP015_PD.global.monmean.prec.nc');
    %ncdisp(path) % check variable name is correct
    preLGM_PI=ncread(path,'aprt'); 

clear path; path=('C:\Users\ph805612\Documents\SISAL\R_programming\Data\Echam5-wiso.T106_EXP016_LGM.global.monmean.tsurf.nc');
    %ncdisp(path) % check variable name is correct
    tmpLGM=ncread(path,'tsurf');
    
clear path; path=('C:\Users\ph805612\Documents\SISAL\R_programming\Data\Echam5-wiso.T106_EXP015_PD.global.monmean.tsurf.nc');
    %ncdisp(path) % check variable name is correct
    tmpLGM_PI=ncread(path,'tsurf');

    
    clear path
    save Runs_LGM
    
    %% 
   clear; close all
   disp('Extracting d18O for LIG runs')

clear path; path=('C:\Users\ph805612\Documents\SISAL\R_programming\Data\Simulating_Climate_and_Stable_Water_Isotopes_during_the_LIG_using_a_Coupled_Climate_Isotope_Model\cosmos-aso-wiso\Eem125-S2\post\echam5\Eem125-S2_echam5_wiso_wisoaprt_d_ymonmean.nc');
    %ncdisp(path) % check variable name is correct
    isoLIG=ncread(path,'wisoaprt_d');
    isoLIG = isoLIG(:,:,2,:); isoLIG = squeeze(isoLIG);
    
    time=ncread(path,'time');
    lon=ncread(path,'lon'); lat=ncread(path,'lat');
    
clear path; path=('C:\Users\ph805612\Documents\SISAL\R_programming\Data\Simulating_Climate_and_Stable_Water_Isotopes_during_the_LIG_using_a_Coupled_Climate_Isotope_Model\cosmos-aso-wiso\EXP003\post\echam5\EXP003_echam5_wiso_wisoaprt_d_ymonmean.nc');
    %ncdisp(path) % check variable name is correct
    isoLIG_PI=ncread(path,'wisoaprt_d'); 
    isoLIG_PI = isoLIG_PI(:,:,2,:); isoLIG_PI=squeeze(isoLIG_PI);

clear path; path=('C:\Users\ph805612\Documents\SISAL\R_programming\Data\Simulating_Climate_and_Stable_Water_Isotopes_during_the_LIG_using_a_Coupled_Climate_Isotope_Model\cosmos-aso-wiso\Eem125-S2\post\echam5\Eem125-S2_echam5_wiso_aprt_ymonmean.nc');
    %ncdisp(path) % check variable name is correct
    preLIG=ncread(path,'aprt'); 
    
clear path; path=('C:\Users\ph805612\Documents\SISAL\R_programming\Data\Simulating_Climate_and_Stable_Water_Isotopes_during_the_LIG_using_a_Coupled_Climate_Isotope_Model\cosmos-aso-wiso\EXP003\post\echam5\EXP003_echam5_wiso_aprt_ymonmean.nc');
    %ncdisp(path) % check variable name is correct
    preLIG_PI=ncread(path,'aprt'); 

clear path; path=('C:\Users\ph805612\Documents\SISAL\R_programming\Data\Simulating_Climate_and_Stable_Water_Isotopes_during_the_LIG_using_a_Coupled_Climate_Isotope_Model\cosmos-aso-wiso\Eem125-S2\post\echam5\Eem125-S2_echam5_main_temp2_ymonmean.nc');
    %ncdisp(path) % check variable name is correct
    tmpLIG=ncread(path,'temp2');    
    
clear path; path=('C:\Users\ph805612\Documents\SISAL\R_programming\Data\Simulating_Climate_and_Stable_Water_Isotopes_during_the_LIG_using_a_Coupled_Climate_Isotope_Model\cosmos-aso-wiso\EXP003\post\echam5\EXP003_echam5_main_temp2_ymonmean.nc');
    %ncdisp(path) % check variable name is correct
    tmpLIG_PI=ncread(path,'temp2'); 
    
time=ncread(path,'time');
lon=ncread(path,'lon'); lat=ncread(path,'lat');

    clear path
    save Runs_LIG
