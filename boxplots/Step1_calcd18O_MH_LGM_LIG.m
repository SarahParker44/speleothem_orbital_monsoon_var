%This script produces the MH.mat, LGM.mat, LIG.mat files with each of the variables: 
%total annual precipitation, mean annual precipitation, temp and
%annual precipitation weighted d18O

clear
close all

for z=1:6
    if z==1

        clearvars -except z
        load Runs_MH
        pre=preMH; iso=isoMH; temp=tmpMH; 
        
    elseif z==2
        
        clearvars -except z
        load Runs_MH
        pre=preMH_PI; iso=isoMH_PI; temp=tmpMH_PI; 
    
    elseif z==3

        clearvars -except z
        load Runs_LGM
        pre=preLGM; iso=isoLGM; temp=tmpLGM; 
        
    elseif z==4
        
        clearvars -except z
        load Runs_LGM
        pre=preLGM_PI; iso=isoLGM_PI; temp=tmpLGM_PI; 
     
    elseif z==5
        
        clearvars -except z
        load Runs_LIG
        pre=preLIG; iso=isoLIG; temp=tmpLIG; 
        
    elseif z==6
        
        clearvars -except z
        load Runs_LIG
        pre=preLIG_PI; iso=isoLIG_PI; temp=tmpLIG_PI;
    end                

        % Calculate weights from field
        pre_res=reshape(pre, size(pre,1)*size(pre,2),size(pre,3));
        evap_res=reshape(evap, size(evap,1)*size(evap,2),size(evap,3));
        temp_res=reshape(temp, size(temp,1)*size(temp,2),size(temp,3));

        weights=ones(size(pre_res,1),size(pre_res,2))*NaN;
        
        pre_annual=ones(size(pre_res,1),size(time,1)/12)*NaN;
        pre_tot_annual=ones(size(pre_res,1),size(time,1)/12)*NaN;
        temp_annual=ones(size(temp_res,1),size(time,1)/12)*NaN;
        
        for st=1:size(pre_res,1) %each gridcell
            if (z==5 || z==6)
                totpre=nansum(pre_res(st,:));
                meanpre=nanmean(pre_res(st,:));
                meantemp=nanmean(temp_res(st,:));
                
                weights(st,:)=pre_res(st,:)/totpre;
                
                pre_annual(st,1)=meanpre;
                pre_tot_annual(st,1)=totpre;
                temp_annual(st,1)=meantemp;
                
            elseif (z==1 || z==2 || z==3 || z==4)
                
                for k=1:size(time,1)/12 %each year
                clear x; x=pre_res(st,12*k-11:12*k); %monthly pre
                clear c; c=temp_res(st,12*k-11:12*k); %monthly temp
                
                %if pre=NaN => convert to 0
                clear m; [m]=find(isnan(x)==1); x(m)=0;
                clear q; [q]=find(isnan(c)==1); c(q)=0;

            totpre=nansum(x); % total mm/month that year
            meanpre=nanmean(x);
            meantemp=nanmean(c);
            
                for j=1:12
                    if x(1,j)==0
                        weights (st,12*k-11+j-1)=0; %ensure that if pre=0 => weight for that month is 0 and not NaN
                    else 
                        weights(st,12*k-11+j-1)=x(1,j)/totpre;
                    end
                end
                
            pre_annual(st,k)=meanpre;
            pre_tot_annual(st,k)=totpre;
            temp_annual(st,k)=meantemp;

                end
            end
        end
        clear st k j totpre totevap y x

        % check that all months with pre=0 have d18O=NaN and not the other
        % way around=> none
        iso_res=reshape(iso,size(iso,1)*size(iso,2),size(iso,3));
        wiso=ones(size(iso_res,1),size(time,1)/12 )*NaN;% Calculate annual weighted d18O means

        for st=1:size(iso_res,1) %each location
            if (z==5 || z==6)
                x=weights(st,:);
                y=iso_res(st,:);
                
                if sum(isnan(y))>0 % this is the number of NaNs in y, if>0, xy=NaN...
                    %remove cols with iso=NaN and then compute wd18O without them
                    clear m; [m]=find(isnan(y)==1);
                    nan_index=fliplr(1:length(m));
                    for j=1:length(nan_index)
                        %remove cols with NaNs, not to be included in the final wd18O
                        %(pre=0 for these cols)
                        y(:,m(nan_index(j)))=[]; x(:,m(nan_index(j)))=[];
                    end
                    xy=sum(x.*y);
                else
                    xy=sum(x.*y); % annual wd18O
                end
                
                wiso(st,1)=xy;
                
            elseif (z==1 || z==2 || z==3 || z==4)
                for k=1:size(time,1)/12 %each year
                    x=weights(st,12*k-11:12*k); % monthly weights
                    y=iso_res(st,12*k-11:12*k); % monthly precip d18O
                    if sum(isnan(y))>0 % this is the number of NaNs in y, if>0, xy=NaN...
                        %remove cols with iso=NaN and then compute wd18O without them
                        clear m; [m]=find(isnan(y)==1);
                        nan_index=fliplr(1:length(m));
                        for j=1:length(nan_index)
                            %remove cols with NaNs, not to be included in the final wd18O
                            %(pre=0 for these cols)
                            y(:,m(nan_index(j)))=[]; x(:,m(nan_index(j)))=[];
                        end
                        xy=sum(x.*y);
                    else
                        xy=sum(x.*y); % annual wd18O
                    end
                    
                    wiso(st,k)=xy;
 
                end
            end
        end

        clear ans st k x y j xy m nan_index 

        pre_ann=reshape(pre_annual,size(lon,1), size(lat,1),size(pre_annual,2));
        pre_tot_ann=reshape(pre_tot_annual,size(lon,1), size(lat,1),size(pre_annual,2));
        temp_ann=reshape(temp_annual,size(lon,1), size(lat,1),size(temp_annual,2));
        wiso=reshape(wiso,size(lon,1), size(lat,1),size(wiso,2));
        pre_mon=reshape(pre_res,size(lon,1), size(lat,1),size(pre_res,2));
        iso_mon=reshape(iso_res,size(lon,1), size(lat,1),size(iso_res,2));
        

        if z==1
            pre_ann_MH=pre_ann; pre_tot_ann_MH=pre_tot_ann; wiso_MH=wiso; temp_MH=temp_ann; 
            clearvars -except wiso_MH temp_MH pre_ann_MH pre_tot_ann_MH mat_MH 
            save('MH.mat')
        elseif z==2
            pre_ann_MH_PI=pre_ann; pre_tot_ann_MH_PI=pre_tot_ann; wiso_MH_PI=wiso;temp_MH_PI=temp_ann;
            clearvars -except wiso_MH_PI temp_MH_PI pre_ann_MH_PI pre_tot_ann_MH_PI 
            save('MH_PI.mat')
        elseif z==3
            pre_ann_LGM=pre_ann; pre_tot_ann_LGM=pre_tot_ann; wiso_LGM=wiso; temp_LGM=temp_ann;
            clearvars -except pre_ann_LGM pre_tot_ann_LGM wiso_LGM temp_LGM 
            save('LGM.mat')
        elseif z==4
            pre_ann_LGM_PI=pre_ann; pre_tot_ann_LGM_PI=pre_tot_ann; wiso_LGM_PI=wiso; temp_LGM_PI=temp_ann;
            clearvars -except pre_ann_LGM_PI pre_tot_ann_LGM_PI wiso_LGM_PI temp_LGM_PI 
            save('LGM_PI.mat')
        elseif z==5
            pre_ann_LIG=pre_ann; pre_tot_ann_LIG=pre_tot_ann; wiso_LIG=wiso; temp_LIG=temp_ann;
            clearvars -except pre_ann_LIG pre_tot_ann_LIG wiso_LIG temp_LIG
            save('LIG.mat')
        elseif z==6
            pre_ann_LIG_PI=pre_ann; pre_tot_ann_LIG_PI=pre_tot_ann; wiso_LIG_PI=wiso; temp_LIG_PI=temp_ann;
            clearvars -except pre_ann_LIG_PI pre_tot_ann_LIG_PI wiso_LIG_PI temp_LIG_PI
            save('LIG_PI.mat')
            
        end  
end
            