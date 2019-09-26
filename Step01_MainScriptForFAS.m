%% Main script for USGS Grant G17AP00034
% Jeff Bayless, AECOM, jeff.bayless@aecom.com
% 26 September 2019

% this script loops through earthquakes within a selected region and does the following:
    % load the data
    % calculate site response factors
    % calculate rad pattern factors
    % adjust the data
    % model the attenuation of the data at 10 frequencies and perform the Q(f) analysis
    % create/save figures and save the results to .mat files

close all; clearvars -except PNUM2 PTXT2 PRAW2; clc;
set(0,'defaultaxesfontsize',12)

%% set some run options

Options.JeffRegionFlag = 2  ; % this is a flag for the region used in selecting the data to analyze.
        % 1 = Gulf Coast, source/site crosses the gulf coast boundary (didn't use this one, used 11 instead)
        % 2 = cena
        % 3 = appalachian province
        % 4 = Atlantic Coastal Plain
        % 5 = all data, irrespective of region
        % 11 = gulf coast; strictly path 1 only
Options.MinMaxDistance = [150 500]; % this is the min and max distance (km) for the data to use in analysis; set to [150 500]
Options.RadType = 2; % pick which radiation pattern adjustment component to use: 1=SH, 2=S (SRSS of SH and SV), 0=off (for testing)
Options.FASAVG = 1; % 1= don't average the EAS, 2= average the EAS in log space; it should be set to 1 since the EAS is already smoothed
Options.CapRad = 0; % 1= cap the radiation pattern adjustment to a factor of 2; 0=no cap
Options.NSigOutlier = 3; % # of sigmas to classify outliers in the attenuation regression, set to a large number to keep all data, using =3
Options.min10 = 10 ; % number of min stats at a freq to use in regressions, set=10
Options.CloseFigs = 1; % set to 1 to close figures after saving them
Options.OtherEASAdj = 1; % set to 1 to calculate c(f) and Q(f) for the EAS without site and radiation pattern adjustments, otherwise these are skipped (leave=1)

Options.SiteAmpModel = 2; % choose which site amplification model from the ones given below

if Options.SiteAmpModel == 1
    SaveDirString='EASRadSite_Harmon';
elseif Options.SiteAmpModel == 2
    SaveDirString='EASRadSite_Stewart';
elseif Options.SiteAmpModel == 3
    SaveDirString='EASRadSite_Bayless';
elseif Options.SiteAmpModel == 4
    SaveDirString='EASRadSite_NoSite';
end
    
mkdir(SaveDirString,['Region' num2str(Options.JeffRegionFlag)]); mkdir([SaveDirString '/Region' num2str(Options.JeffRegionFlag) '/Figures']);

%% load some files

% state and lake outlines
Mstate=m_shaperead('/Users/jeff/WORK/MATLAB/m_map/ne_50m_admin_1_states_provinces_lakes/ne_50m_admin_1_states_provinces_lakes');
Mlake=m_shaperead('/Users/jeff/WORK/MATLAB/m_map/ne_50m_lakes/ne_50m_lakes');
% PEER 2014 boundaries
[PEERceus,~,~]=xlsread('AppendixE-PEER-2014-15.xlsx','Sheet1','A3:B239'); % Central US, Canada (CNA)	
[PEERgulf,~,~]=xlsread('AppendixE-PEER-2014-15.xlsx','Sheet1','D3:E197'); % Gulf Coast and Mississippi Embayment (MEM)	
[PEERapp,~,~]=xlsread('AppendixE-PEER-2014-15.xlsx','Sheet1','G3:H747'); % Appalachians (APP)	
[PEERacp,~,~]=xlsread('AppendixE-PEER-2014-15.xlsx','Sheet1','J3:K674'); % Atlantic Coastal Plain (ACP)	
% Cramer (2014) boundaries
load CramerBounds.mat %,'Cram1','Cram1b','Cram2')
% PEER (2015) Chapter 3 crustal model
%   thickness depth   Vs (km/s)
PCrust=[1.0	 1.0 3.00
        11.0 12.0	3.52
        28.0 40.0	3.80];

%% read the PEER flatfile, load EAS
if exist('PNUM2','var')
    % continue
else
    [PNUM2,PTXT2,PRAW2] = xlsread('/Users/jeff/WORK/USGS/USGS_2017/CEUS Q/Goulet_Bayless/NGA-East_RotD50_5pct_Flatfile_Public_20141118_FASadded.xlsx','NGA-East_RotD50_5pct_v20141118','A2:XK9383');  
end

EQName=PTXT2(:,1); M=PNUM2(:,10);
RSN=PNUM2(:,1); EQID=PNUM2(:,2);
EpiLat=PNUM2(:,12); EpiLon=PNUM2(:,13);
HypDep=PNUM2(:,16); Strike=PNUM2(:,17);
Dip=PNUM2(:,18); Rake=PNUM2(:,19);
Mech=PNUM2(:,21); PIE=PNUM2(:,22);
SSN=PNUM2(:,13); StaName=PRAW2(:,29);
SLat=PNUM2(:,30); SLon=PNUM2(:,31);
Vs30=PNUM2(:,34); Repi=PNUM2(:,37);
Rhyp=PNUM2(:,38); Rrup=PNUM2(:,39);
Rjb=PNUM2(:,40); PGA=PNUM2(:,86);
EvRegion=PNUM2(:,41); StaRegion=PNUM2(:,42); PathRegion=PNUM2(:,43);

iFreqs=[518 519   534   546   559   571   584   597   611   623   ];  % pick the columns of FAS data at the freqs given by FREQS
FREQS=[1.0 1.51 2.14 2.81 3.80 5.01 6.76 9.12 12.59 16.6 ];

if Options.FASAVG==1 % use the EAS(f) directly
    FAS=PNUM2(:,iFreqs);
    ind99=FAS<0;
    FAS(ind99)=nan;
elseif Options.FASAVG==2    % take the average FAS over a small frequency band
    for ii=1:length(iFreqs)
        if ii<3
            RNG=PNUM2(:,iFreqs(ii):iFreqs(ii)+7);
        else
            RNG=PNUM2(:,iFreqs(ii)-7:iFreqs(ii)+7);
        end
        indrng=RNG<0;
        RNG(indrng)=nan;
        FAS(:,ii)=exp(nanmean(log(RNG),2));
    end
    infasp=FAS<0; FAS(infasp)=nan;
end

%% Event selection
if Options.JeffRegionFlag==11 % Gulf Coast only (strict)
    GCEQID=[46	47	48	49	66	67	73	74	76	80	81	90	91	92 760 7600 800 8000]; % didn't use: 34 38 41 58
    EvUnique=GCEQID;
elseif Options.JeffRegionFlag==2 % CENA
    CENAID=[16 25 29 30 33 44 46 47  55 56 57 58 60 61 66 67 73 74 75 76 80 81 83 90 91]; % didnt use: 21 32 34 35 37 51 85 48
    EvUnique=CENAID;
elseif Options.JeffRegionFlag==3 % Appalachian
    APPID=[16 26 51 68 82 85 88 89 ]; % didnt use: 40 60 61 64 21
    EvUnique=APPID;
elseif Options.JeffRegionFlag==4 % atlantic coast   
    EvUnique=[21 88]; % skip 20 80
    Options.min10=3;
    Options.MinMaxDistance = [100 600];
elseif Options.JeffRegionFlag==5 % all data
    [EvUnique,IA,IC] = unique(EQID,'stable');
    % keep only the events that have enough data, etc. these are hand-picked
    keepthese=[4 5 10 11 12 14 16 17 18 19 21 22 24 25 26 27 28 29 30 32 33 34 35 36 38 39 40 41 43 47 48 49 50 52 53 55 57 58 59 63 64 65 66 67 68 69 70 72 75 76 77 78 79];
    EvUnique=EvUnique(keepthese);
    IA=IA(keepthese);
end

%% Main loop
Qall=nan(length(EvUnique),10);
QallPlus=nan(length(EvUnique),10);
QallMinus=nan(length(EvUnique),10);

for ii=1:length(EvUnique) 
    
    Evi=EvUnique(ii);
    
    if Evi==760 % Guy eqk, east
        SpecialEventFlag=1;
        Evi=76;
    elseif Evi==7600 % Guy eqk, west
        SpecialEventFlag=2;
        Evi=76;
    elseif Evi==800 % Greenbrier eqk, east
        SpecialEventFlag=3;
        Evi=80;
    elseif Evi==8000 % Greenbier eqk, west
        SpecialEventFlag=4;
        Evi=80;
    elseif Evi==47 && Options.JeffRegionFlag==2 % Mt carmel
        SpecialEventFlag=11;
    elseif Evi==55 && Options.JeffRegionFlag==2 % constance bay, no strike/dip info so no radpat adjust
        SpecialEventFlag=12;
    elseif Evi==83 && Options.JeffRegionFlag==2 % val de bois, no strike/dip info so no radpat adjust
        SpecialEventFlag=13;
    elseif Evi==25 && Options.JeffRegionFlag==2 % bark lake, no strike/dip info so no radpat adjust
        SpecialEventFlag=13;
    elseif Evi==51 && Options.JeffRegionFlag==3 % Riviere du loup; S radiation pattern fits better
        SpecialEventFlag=6;
    elseif Evi==89 && Options.JeffRegionFlag==3 % Mineral; S radiation pattern fits better
        SpecialEventFlag=6;
    elseif Evi==16 && Options.JeffRegionFlag==3 % Aus Sable; modify Qpos selection
        SpecialEventFlag=7;
    elseif Evi==26 && Options.JeffRegionFlag==3 % Jefferson; modify data selection (path 3 only)
        SpecialEventFlag=8;
    elseif Evi==82 && Options.JeffRegionFlag==3 % Eagle lake; same as Jefferson
        SpecialEventFlag=8;
    elseif Evi==82 && Options.JeffRegionFlag==3 % Hawkesbury; extend to 600 km to get data
        SpecialEventFlag=9;
    elseif Evi==68 && Options.JeffRegionFlag==3 % concord, no strike/dip info so no radpat adjust
        SpecialEventFlag=13;
    elseif Evi==88 && Options.JeffRegionFlag==4 % Mineral, Qposi reqts
        SpecialEventFlag=10;
    else
        SpecialEventFlag=0;
    end
    EvInd=EQID==Evi;
    nReci=sum(EvInd);
    i1=find(EvInd==1,1);% index of the first row for this event in the flatfile
    EvNamei=EQName{i1};
    
    % using indexing, keep the recs for this event, after applying filters
    indStaReg1= StaRegion==1 | PathRegion==1 | PathRegion==6; % gulf coast; PathRegion6 means source/site crosses the gulf coast boundary
    indStaReg1strict= StaRegion==1 | PathRegion==1 | PathRegion==1; % gulf coast; path 1 only
    indStaReg2= StaRegion==2 | PathRegion==2 ;%| PathRegion==5; % cena
    indStaReg3= StaRegion==3 | PathRegion==3 | PathRegion==5; % appalachian province
    indStaReg4= StaRegion==4 | PathRegion==4 ;%| PathRegion==5; % Atlantic Coastal Plain
    indStaRegAll= StaRegion>0 & PathRegion>0;

    if Options.JeffRegionFlag==1        
        % 1 = Gulf Coast, source/site crosses the gulf coast boundary
        indStaRegUSE=indStaReg1;
        RegionString='Gulf Coast';
    elseif Options.JeffRegionFlag==2   
        % 2 = cena
        indStaRegUSE=indStaReg2;
        RegionString='CENA';
    elseif Options.JeffRegionFlag==3  
        % 3 = appalachian province
        indStaRegUSE=indStaReg3;
        RegionString='Appalachian Province';
    elseif Options.JeffRegionFlag==4   
        % 4 = Atlantic Coastal Plain
        indStaRegUSE=indStaReg4;
        RegionString='Atlantic Coastal Plain';
    elseif Options.JeffRegionFlag==5   
        % 5 = all data, irrespective of region
        indStaRegUSE=indStaRegAll;
        RegionString='All';
    elseif Options.JeffRegionFlag==11   
        % 11 = gulf coast; strictly path 1 only (no crossing into)
        indStaRegUSE=indStaReg1strict;
        RegionString='Gulf Coast, Strict';
    end
    RegString=['_Region' num2str(Options.JeffRegionFlag)];
    
    Rmax=Options.MinMaxDistance(2); % max distance for data
    Rmin=Options.MinMaxDistance(1); % min
    if SpecialEventFlag==8
        Rmax=600;
        indStaRegUSE=StaRegion==3 | PathRegion==3 ;
    end
    if SpecialEventFlag==9
        Rmax=600;
    end
    if SpecialEventFlag==11
        Rmax=450;
        Options.min10=9;
    end

    indR = Rrup>Rmin & Rrup<Rmax;
    ind = EvInd & indStaRegUSE & indR ;

    if sum(ind)<Options.min10
        disp(['Event ' num2str(Evi) ', ' EvNamei ',  has only ' num2str(sum(ind),'%d') ' out of ' num2str(nReci) ' compatible stats, skipping'])
        continue
    end

    Strikei0=Strike(ind); Rakei0=Rake(ind);
    Dipi0=Dip(ind); Vi0=Vs30(ind);
    Ri0=Rrup(ind); Mi0=M(ind);
    HypD0=HypDep(ind);
    RSNi=RSN(ind); FASi=FAS(ind,:);
    StaLL=[SLon(ind) SLat(ind)];
    
    % setup map
    f1=figure(1); set(f1,'position',[1   156   832   652]);%[1         156        1031         652]);
    m_proj('lambert','long',[-108 -60],'lat',[24 56]); 
    [CS,CH]=m_etopo2('contourf',[-10000 3:250:4000],'edgecolor','none');
    colormap([repmat(rgb('royalblue'),125,1); m_colmap('gland',80)]);
    brighten(.5);
    % state outlines
    for k=1:size(Mstate.ncst,1)
         m_plot(Mstate.ncst{k}(:,1),Mstate.ncst{k}(:,2),'color',rgb('grey'),'linewidth',.025); 
    end
    % lakes
    for k=[5 10 12 16 22 23 24 25]
            if k==25
                mind=1:size(Mlake.ncst{k},1)-11;
            else
                mind=1:size(Mlake.ncst{k},1);
            end
            ml=m_patch(Mlake.ncst{k}(mind,1),Mlake.ncst{k}(mind,2),[.4 .58 .93],'edgecolor','none'); %rgb('cornflowerblue')        
    end
    % plot the epicenter
    m_plot(EpiLon(i1),EpiLat(i1),'kp','markerfacecolor','r','markersize',9); hold on

    SH=zeros(size(FASi,2),360);
    fhdist=figure(100+ii); set(fhdist,'position',[360         172        1035         742]);%[1    66   988   919])
    ha = tight_subplot(3,3,[.02 .02],[.05],[.05]);
    for f=1:size(FASi,2) % loop over the set of ~10 frequencies

        indFAS=isfinite(FASi(:,f));
        
        % for events like Guy, with clear breaks in data regions, investigate differences
        if SpecialEventFlag==1
            indFAS=isfinite(FASi(:,f)) & StaLL(:,1)>-91.1 & StaLL(:,2)>33.4;
            Evi=760;
        elseif SpecialEventFlag==2
            indFAS=isfinite(FASi(:,f)) & StaLL(:,2)<=33.4;
            Evi=7600;
        elseif SpecialEventFlag==3
            indFAS=isfinite(FASi(:,f)) & StaLL(:,1)>-91.1 & StaLL(:,2)>33.4;
            Evi=800;
        elseif SpecialEventFlag==4
            indFAS=isfinite(FASi(:,f)) & StaLL(:,2)<=33.4;
            Evi=8000;
        end
        if sum(indFAS)<Options.min10
            disp(['Event ' num2str(Evi) ', ' EvNamei ',  has only ' num2str(sum(indFAS),'%d') ' compatible stats, skipping f = ' num2str(FREQS(f))])
            continue
        end

        Ri=Ri0(indFAS); Vi=Vi0(indFAS);
        Mi=Mi0(1); Strikei=Strikei0(1);
        Dipi=Dipi0(1); Rakei=Rakei0(1);
        HD=HypD0(1);
        LLi=StaLL(indFAS,:);
        fi=FREQS(f); % Hz  
        if fi==5.01
            Ritrack=Ri;
        end
        Yraw=log(FASi(indFAS,f)); % the ground motion at this freq, ln units

        %% adjustment for site            
        % use: Harmon, et al. 2019
        AMP_Harmon=zeros(size(Vi));
        for vv=1:length(Vi)
            AMP_Harmon(vv,1) = Harmon2019SiteAmpLinear(fi,Vi(vv));
        end
        % Stewart et al. (2017)  (for PSA)
        FvT=zeros(size(Vi));
        for vv=1:length(Vi)
            [~,FvT(vv,1)]=ceus_site_lin_Fv(Vi(vv),1./fi);
        end
        Ysitetemp=Yraw-FvT;
        % Bayless and Abrahamson (2019) (for WUS)
        FBA=zeros(size(Vi));
        for vv=1:length(Vi)
            [FBA(vv,1)]=Bayless2019SiteAmp(fi,Vi(vv),760);
        end

        %% adjustment for radiation pattern
        if fi<3 %the rad pattern washes out at f>3 Hz

            f1rad=1; f2rad=3; 
            [SH(f,:),SV(f,:),~,S(f,:)]=rad_pat(Strikei,Dipi,Rakei,f1rad,f2rad,fi); % get the theoretical rad pattern at this freq for tko=120
            if Options.RadType==1
                RGeo=SH(f,:);
            elseif Options.RadType==2
                RGeo=S(f,:);
            elseif Options.RadType==0
                RGeo=ones(size(SH(f,:)));
            end
            if SpecialEventFlag==12 || SpecialEventFlag==13
                RGeo=ones(size(SH(f,:)));
            end
            if f==2  % plot the radiation pattern at 1.5 Hz
                figure(1)
                [e1,e2]=m_ll2xy(EpiLon(i1),EpiLat(i1));
                azr2=1:360;
                figure(1)
                prad=plot(RGeo/10.*cosd(azr2)+e1,RGeo/10.*sind(azr2)+e2,'--','color',rgb('salmon'),'linewidth',1.5);          
                m_scatter(LLi(:,1),LLi(:,2),50./max(exp(Ysitetemp)).*exp(Ysitetemp),'k^','filled','markerfacecolor',rgb('teal'),'markeredgecolor','k'); hold on 
            end
            if SpecialEventFlag==6 || SpecialEventFlag==8 || SpecialEventFlag==11
                RGeo=S(f,:);
            end
            if SpecialEventFlag==13
                RGeo=SH(f,:);
            end

            distcheck=zeros(size(Vi)); brng=zeros(size(Vi)); Yrad=zeros(size(Vi)); Radtrack=zeros(size(Vi)); aztrack=zeros(size(Vi)); RadLog=zeros(size(Vi)); aztrack2=zeros(size(Vi));
            for ss=1:length(Vi) % loop over the sites
                [distcheck(ss,1),brng(ss,1),~,~] = calc_dist_brng(EpiLat(i1),EpiLon(i1),LLi(ss,2),LLi(ss,1));
                azEN=brng(ss,1); % this is source-site azimuth in degrees east of north   
                azindex=ceil(wrapTo360(90-azEN)); % this is the angle measured CCW from East
                % Note for different notations in Matlab and e.g. Aki&Richards:
                % Matlab definition of spherical coordinates [R,THETA,PHI]:
                %  "THETHA is the counterclockwise angle in the xy plane measured from the
                %    positive x axis.  PHI is the elevation angle from the xy plane.
                %
                % Common (e.g. Aki&Richards) definition:
                %  "phi is the clockwise angle in the xy plane measured from the
                %    positive x axis.  theta is the vertical angle starting from
                %    the y axis"

                RadAdjustment=RGeo(azindex)/mean(RGeo); % this is the ratio of the rad pattern at this azimuth to the average over all azimuths
                if Options.CapRad==1
                    if RadAdjustment<0.5
                        disp([' Station number ' num2str(ss) ' has RadPattern Adjustment ' num2str(RadAdjustment) ])
                        RadAdjustment=0.5;
                    elseif RadAdjustment>2
                        disp([' Station number ' num2str(ss) ' has RadPattern Adjustment ' num2str(RadAdjustment) ])
                        RadAdjustment=2;
                    end
                end
                Radtrack(ss,1)=RadAdjustment;
                aztrack(ss,1)=azindex;                
                aztrack2(ss,1)=ceil(wrapTo360(ceil(azEN)-Strikei));
                RadLog(ss,1)=log(RadAdjustment); % ln units
            end

        else
            RadLog=zeros(size(Vi)); % no rad pattern adjustment for freq>=3 Hz
        end
        if f==2
            % plot the site adjustment
            fsit=figure(800+ii); set(fsit,'position',[360   595   601   319]);%[1    66   988   919])
            loglog(Vi,exp(FvT),'ko','markerfacecolor','r'); hold on
            plot([200 1500],[1 1],'k-','linewidth',.7)
            grid on; 
            axis([200 1500 0.4 2])
            xlabel('V_{s30} (m/s)')
            ylabel('Sea17 Site Response Adjustment, \itF_{Site}')

            xl = get(gca,'XLabel'); yl = get(gca,'YLabel');
            xlFontSize = get(xl,'FontSize');
            xAX = get(gca,'XAxis');
            set(xAX,'FontSize', 9)
            yAX = get(gca,'YAxis');
            set(yAX,'FontSize', 6);
            set(xl, 'FontSize', xlFontSize);
            set(yl, 'FontSize', xlFontSize);

            text(220,0.5,[' f = ' num2str(fi,'%.2f') ' Hz'],'edgecolor','k','backgroundcolor','w')
            jeffsavefig_r300(fsit,[SaveDirString '/Region' num2str(Options.JeffRegionFlag) '/Figures/' num2str(Evi,'%d') '_SiteAdj' RegString '.png']);
            if Options.CloseFigs
                close(fsit)
            end
        end

        %% the corrected data to use in regression
        
        Yrad=Yraw-RadLog; % this is EAS_rad   
        switch SaveDirString
            case 'EASRadSite_Stewart'
                Ysite=Yraw-FvT; % this is EAS_site using the Stewart correction
                Yradsite=Yraw-FvT-RadLog; % this is EAS_radsite 
            case 'EASRadSite_Bayless'
                Ysite=Yraw-FBA; % this is EAS_site using the Bayless correction
                Yradsite=Yraw-FBA-RadLog; % this is EAS_radsite 
            case 'EASRadSite_Harmon'
                Ysite=Yraw-AMP_Harmon; % this is EAS_site  using Harmon  
                Yradsite=Yraw-AMP_Harmon-RadLog; % this is EAS_radsite 
            case 'EASRadSite_NoSite'
                Ysite=Yraw; %  
                Yradsite=Yrad; % this is EAS_radsite 
        end
        % pick the one for plotting the results
        Y=Yradsite;     
        
        %% model the attenuation with distance

        options = statset('TolX',1e-15,'TolFun',1e-15,'robust','on');
        options.RobustWgtFun = 'huber'; %'bisquare';
        beta0=[max(Y);-.005]; % initial values for the fit
        [beta1,betares,jaco,covb,mse]=nlinfit(Ri,Y,@FASFIT,beta0,options); % nonlinear fit using Cramer (2014) eqn 1
        xdo=min(Ri):5:max(Ri); % range of R values
        [ydo,ydelta] = nlpredci(@FASFIT,xdo.',beta1,betares,'covar',covb); % prediction for these R, with 95% CI
        ysig=std(betares); % sigma of the residuals
        YSIG(ii,f)=ysig; betares0=betares;
        
        % track and exclude any outliers (more than Options.NSigOutlier sigma away)
        nsigoutlier=Options.NSigOutlier;
        outliers = abs(betares)>nsigoutlier*ysig;
        noout= abs(betares)<=nsigoutlier*ysig;
        if sum(outliers)>0
            disp([ num2str(sum(outliers)) ' outliers!!!'])
            beta0=[max(Y);-.015]; % initial values for the fit
            [beta1,betares,jaco,covb,mse]=nlinfit(Ri(noout),Y(noout),@FASFIT,beta0,options); % nonlinear fit using Cramer (2014) eqn 1
            [ydo,ydelta] = nlpredci(@FASFIT,xdo.',beta1,betares,'covar',covb); % prediction for these R, with 95% CI
            ysig=std(betares); % sigma of the residuals
            YSIG(ii,f)=ysig;
        end

        %% calculate Q
        Cevent=beta1(2); % C from Cramer (2014) eqn 1
        CeventSigma=sqrt(covb(2,2));
%         Beta=3.5; %km/s shear wave velocity near source
        Beta=interp1(PCrust(:,2),PCrust(:,3),HD); 
        Qevent=-pi.*fi./Cevent./Beta; % apparent Q, Cramer (2014) eqn 2
        Call(ii,f)=Cevent;
        Qall(ii,f)=Qevent;
        % plus or minus one sigma
        QallPlus(ii,f)=Qevent.*Cevent./(Cevent+CeventSigma); % this is -pi*f/beta/(c+sigma)
        QallMinus(ii,f)=Qevent.*Cevent./(Cevent-CeventSigma);
        
        CallPlus(ii,f)=(Cevent+CeventSigma);
        CallMinus(ii,f)=(Cevent-CeventSigma);
        
        if Options.OtherEASAdj % do the fit for EAS_raw and EAS_Site
            
            [betaS,betaresS,jacoS,covbS,mseS]=nlinfit(Ri,Ysite,@FASFIT,beta0,options); % fitting the EAS_site data
            [betaR,betaresR,jacoR,covbR,mseR]=nlinfit(Ri,Yrad,@FASFIT,beta0,options); % fitting EAS_rad
            [betaRaw,betaresRaw,jacoRaw,covbRaw,mseRaw]=nlinfit(Ri,Yraw,@FASFIT,beta0,options); % fitting EAS_raw
            YSIGS(ii,f)=std(betaresS); % sigma of the residuals
            YSIGR(ii,f)=std(betaresR); % sigma of the residuals
            YSIGRAW(ii,f)=std(betaresRaw); % sigma of the residuals
            
            % calculate Q for EAS_site
            CeventS=betaS(2); % C from Cramer (2014) eqn 1
            CeventSigmaS=sqrt(covbS(2,2));
            
            QeventS=-pi.*fi./CeventS./Beta; % apparent Q, Cramer (2014) eqn 2
            CallS(ii,f)=CeventS;
            QallS(ii,f)=QeventS;
            % plus or minus one sigma
            QallPlusS(ii,f)=QeventS.*CeventS./(CeventS+CeventSigmaS); % this is -pi*f/beta/(c+sigma)
            QallMinusS(ii,f)=QeventS.*CeventS./(CeventS-CeventSigmaS);
            CallPlusS(ii,f)=(CeventS+CeventSigmaS);
            CallMinusS(ii,f)=(CeventS-CeventSigmaS);

            
            % calculate Q for EAS_rad
            CeventR=betaR(2); % C from Cramer (2014) eqn 1
            CeventSigmaR=sqrt(covbR(2,2));
            
            QeventR=-pi.*fi./CeventR./Beta; % apparent Q, Cramer (2014) eqn 2
            CallR(ii,f)=CeventR;
            QallR(ii,f)=QeventR;
            % plus or minus one sigma
            QallPlusR(ii,f)=QeventR.*CeventR./(CeventR+CeventSigmaR); % this is -pi*f/beta/(c+sigma)
            QallMinusR(ii,f)=QeventR.*CeventR./(CeventR-CeventSigmaR);
            CallPlusR(ii,f)=(CeventR+CeventSigmaR);
            CallMinusR(ii,f)=(CeventR-CeventSigmaR);
            
            % calculate Q for EAS_raw
            CeventRaw=betaRaw(2); % C from Cramer (2014) eqn 1
            CeventSigmaRaw=sqrt(covbRaw(2,2));
            
            QeventRaw=-pi.*fi./CeventRaw./Beta; % apparent Q, Cramer (2014) eqn 2
            CallRaw(ii,f)=CeventRaw;
            QallRaw(ii,f)=QeventRaw;
            % plus or minus one sigma
            QallPlusRaw(ii,f)=QeventRaw.*CeventRaw./(CeventRaw+CeventSigmaRaw); % this is -pi*f/beta/(c+sigma)
            QallMinusRaw(ii,f)=QeventRaw.*CeventRaw./(CeventRaw-CeventSigmaRaw);
            CallPlusRaw(ii,f)=(CeventRaw+CeventSigmaRaw);
            CallMinusRaw(ii,f)=(CeventRaw-CeventSigmaRaw);
        end

        %% distance scaling plot
        if f>1
            figure(100+ii);
            axes(ha(f-1)); 
            loglog(Ri(noout),exp(Y(noout)),'ko','markerfacecolor',rgb('royalblue')); hold on
            if sum(outliers)>0
                loglog(Ri(outliers),exp(Y(outliers)),'kx','markerfacecolor',rgb('black'));
            end
            grid on; hold on
            loglog(xdo,exp(ydo),'r-','linewidth',2)
            loglog(xdo,exp(ydo+ysig),'r--','linewidth',2)
            loglog(xdo,exp(ydo-ysig),'r--','linewidth',2)
            loglog(xdo,exp(ydo(1)).*sqrt(xdo(1))./sqrt(xdo),'k-','linewidth',1.5)
            if f>7
            xlabel('Rupture Distance (km)')
            end
            if f==2 || f==5 || f==8
            ylabel(['EAS (g-s)'])
            end
            yma=ceil(log10(max(exp(Y))));
            ymi=yma-3; 
            axis([100 600 10^ymi 10^yma])
              
            if f==3
            title(['M' num2str(Mi,'%.1f') ' ' strrep(EvNamei,'_',' ') ', EQID ' num2str(Evi,'%d') ', Region: ' RegionString ])
            end
            if f==2
                text(75,1.8.*10^yma,'(b)','fontsize',20);
            end
            tstr=['\itc\rm = ' num2str(Cevent,'%.4f') ', \itQ\rm = ' num2str(round(Qevent),'%d')];
            yt=2*10^ymi;

            xl = get(gca,'XLabel'); yl = get(gca,'YLabel');
            xlFontSize = get(xl,'FontSize');
            xAX = get(gca,'XAxis');
            set(xAX,'FontSize', 9)
            yAX = get(gca,'YAxis');
            set(yAX,'FontSize', 6);
            set(xl, 'FontSize', xlFontSize);
            set(yl, 'FontSize', xlFontSize);

            text(110,yt,tstr,'edgecolor','k','backgroundcolor','w')
            text(350,yt,[' f = ' num2str(fi,'%.2f') ' Hz'],'edgecolor','k','backgroundcolor','w')
        end
        
        if f==3
            %% distance, az, Vs30 resids plot
            fresid=figure(500+ii); set(fresid,'position',[360         172        1035         742]);%[1    66   988   919])
            ha2 = tight_subplot(3,3,[.02 .02],[.05],[.05]);
            figure(500+ii);
                axes(ha2(1)); % distance dep of resids before correcting
                semilogx(Ri,betaresRaw,'ko','markerfacecolor',rgb('royalblue'),'markersize',11); hold on
                grid on; hold on
                ylabel('ln residual, \itEAS\rm_{raw}')
                axis([100 600 -2 2])
                xl = get(gca,'XLabel'); yl = get(gca,'YLabel');
                xlFontSize = get(xl,'FontSize');
                xAX = get(gca,'XAxis');
                set(xAX,'FontSize', 11)
                yAX = get(gca,'YAxis');
                set(yAX,'FontSize', 11);
                set(xl, 'FontSize', xlFontSize);
                set(yl, 'FontSize', xlFontSize);
    %             text(85,1.9,'(a)','fontsize',18);
                tstr=char(['\sigma_{Resid} = ' num2str(std(betaresRaw),'%.2f') ],['se_c = ' num2str(CeventSigmaRaw,'%.2e')]);
                text(105,1.4,tstr,'edgecolor','k','backgroundcolor','w','fontsize',12)

            axes(ha2(4)); % distance depndence of resids after correction
                semilogx(Ri,betaresS,'ko','markerfacecolor',rgb('indianred'),'markersize',11); hold on
                grid on; hold on
                ylabel('ln residual, \itEAS\rm_{Site}')
                axis([100 600 -2 2])
    %             text(85,1.9,'(c)','fontsize',18);
                tstr=char(['\sigma_{Resid} = ' num2str(std(betaresS),'%.2f') ],['se_c = ' num2str(CeventSigmaS,'%.2e')]);
                text(105,1.4,tstr,'edgecolor','k','backgroundcolor','w','fontsize',12)
                xl = get(gca,'XLabel'); yl = get(gca,'YLabel');
                xlFontSize = get(xl,'FontSize');
                xAX = get(gca,'XAxis');
                set(xAX,'FontSize', 10)
                yAX = get(gca,'YAxis');
                set(yAX,'FontSize', 11);
                set(xl, 'FontSize', xlFontSize);
                set(yl, 'FontSize', xlFontSize);

            axes(ha2(7)); % distance depndence of resids after correction
                semilogx(Ri,betares0,'ko','markerfacecolor',rgb('green'),'markersize',11); hold on
                grid on; hold on
                xlabel('Rupture Distance (km)')
                ylabel('ln residual, \itEAS\rm_{RadSite}')
                axis([100 600 -2 2])
    %             text(85,1.9,'(c)','fontsize',18);
                tstr=char(['\sigma_{Resid} = ' num2str(std(betares0),'%.2f') ],['se_c = ' num2str(CeventSigma,'%.2e')]);
                text(105,1.4,tstr,'edgecolor','k','backgroundcolor','w','fontsize',12)
                xl = get(gca,'XLabel'); yl = get(gca,'YLabel');
                xlFontSize = get(xl,'FontSize');
                xAX = get(gca,'XAxis');
                set(xAX,'FontSize', 10)
                yAX = get(gca,'YAxis');
                set(yAX,'FontSize', 11);
                set(xl, 'FontSize', xlFontSize);
                set(yl, 'FontSize', xlFontSize);

                % azimuth

            axes(ha2(2)); % azimuthal dep of resids before correcting
                plot(aztrack2,betaresRaw,'kd','markerfacecolor',rgb('royalblue'),'markersize',11); hold on
                grid on; hold on
                axis([0 360 -2 2])
                xl = get(gca,'XLabel'); yl = get(gca,'YLabel');
                xlFontSize = get(xl,'FontSize');
                xAX = get(gca,'XAxis');
                set(xAX,'FontSize', 11)
                yAX = get(gca,'YAxis');
                set(yAX,'FontSize', 11);
                set(xl, 'FontSize', xlFontSize);
                set(yl, 'FontSize', xlFontSize);
    %             text(-16,1.9,'(b)','fontsize',18);

            axes(ha2(5)); % distance depndence of resids after correction
                plot(aztrack2,betaresS,'kd','markerfacecolor',rgb('indianred'),'markersize',11); hold on
                grid on; hold on
                axis([0 360 -2 2])
                xl = get(gca,'XLabel'); yl = get(gca,'YLabel');
                xlFontSize = get(xl,'FontSize');
                xAX = get(gca,'XAxis');
                set(xAX,'FontSize', 10)
                yAX = get(gca,'YAxis');
                set(yAX,'FontSize', 11);
                set(xl, 'FontSize', xlFontSize);
                set(yl, 'FontSize', xlFontSize);
    %             text(-16,1.9,'(d)','fontsize',18);

            axes(ha2(8)); % distance depndence of resids after correction
                plot(aztrack2,betares0,'kd','markerfacecolor',rgb('green'),'markersize',11); hold on
                grid on; hold on
                xlabel('Source-Site Azimuth (deg)')
                axis([0 360 -2 2])
                xl = get(gca,'XLabel'); yl = get(gca,'YLabel');
                xlFontSize = get(xl,'FontSize');
                xAX = get(gca,'XAxis');
                set(xAX,'FontSize', 10)
                yAX = get(gca,'YAxis');
                set(yAX,'FontSize', 11);
                set(xl, 'FontSize', xlFontSize);
                set(yl, 'FontSize', xlFontSize);
    %             text(-16,1.9,'(d)','fontsize',18);

                % Vs30

            axes(ha2(3)); % azimuthal dep of resids before correcting
                semilogx(Vi,betaresRaw,'k^','markerfacecolor',rgb('royalblue'),'markersize',11); hold on
                grid on; hold on
                axis([200 1500 -2 2])
                xl = get(gca,'XLabel'); yl = get(gca,'YLabel');
                xlFontSize = get(xl,'FontSize');
                xAX = get(gca,'XAxis');
                set(xAX,'FontSize', 11)
                yAX = get(gca,'YAxis');
                set(yAX,'FontSize', 11);
                set(xl, 'FontSize', xlFontSize);
                set(yl, 'FontSize', xlFontSize);
    %             text(60,1.9,'(b)','fontsize',18);

            axes(ha2(6)); % distance depndence of resids after correction
                semilogx(Vi,betaresS,'k^','markerfacecolor',rgb('indianred'),'markersize',11); hold on
                grid on; hold on
                axis([200 1500 -2 2])
                xl = get(gca,'XLabel'); yl = get(gca,'YLabel');
                xlFontSize = get(xl,'FontSize');
                xAX = get(gca,'XAxis');
                set(xAX,'FontSize', 10)
                yAX = get(gca,'YAxis');
                set(yAX,'FontSize', 11);
                set(xl, 'FontSize', xlFontSize);
                set(yl, 'FontSize', xlFontSize);
    %             text(60,1.9,'(d)','fontsize',18);

            axes(ha2(9)); % distance depndence of resids after correction
                semilogx(Vi,betares0,'k^','markerfacecolor',rgb('green'),'markersize',11); hold on
                grid on; hold on
                xlabel('V_{s30} (m/s)')
                axis([200 1500 -2 2])
                xl = get(gca,'XLabel'); yl = get(gca,'YLabel');
                xlFontSize = get(xl,'FontSize');
                xAX = get(gca,'XAxis');
                set(xAX,'FontSize', 10)
                yAX = get(gca,'YAxis');
                set(yAX,'FontSize', 11);
                set(xl, 'FontSize', xlFontSize);
                set(yl, 'FontSize', xlFontSize);
    %             text(60,1.9,'(d)','fontsize',18);

            set(ha2(1:6),'XTickLabel',''); %set(ha(hai),'YTickLabel','')
            set(ha2([2 3  5 6 8 9]),'YTickLabel',''); %set(ha(hai),'YTickLabel','')
            set(ha2,'XMinorGrid','off')
            set(ha2,'XMinorTick','off')

            jeffsavefig_r300(fresid,[SaveDirString '/Region' num2str(Options.JeffRegionFlag) '/Figures/' num2str(Evi,'%d') '_Resids' RegString '.png']);
            if Options.CloseFigs
                close(fresid)
            end
            
            if Evi==90
%        plot attenuation at this freq alone for the example in the report
                fat2=figure(600+ii); set(fat2,'position',[360   595   601   319]);%[1    66   988   919])
                loglog(Ri,exp(Y),'ko','markerfacecolor',rgb('royalblue')); hold on
                grid on; hold on
                loglog(xdo,exp(ydo),'r-','linewidth',2)
                loglog(xdo,exp(ydo+ysig),'r--','linewidth',2)
                loglog(xdo,exp(ydo-ysig),'r--','linewidth',2)
                loglog(xdo,exp(ydo(1)).*sqrt(xdo(1))./sqrt(xdo),'k-','linewidth',1.5)
                xlabel('Rupture Distance (km)')
                ylabel(['\itEAS_{RadSite}\rm (g-s)'])
                yma=ceil(log10(max(exp(Y))));
                ymi=yma-2; 
                axis([100 600 10^ymi 10^yma])
                yt=2*10^ymi;

                xl = get(gca,'XLabel'); yl = get(gca,'YLabel');
                xlFontSize = get(xl,'FontSize');
                xAX = get(gca,'XAxis');
                set(xAX,'FontSize', 9)
                yAX = get(gca,'YAxis');
                set(yAX,'FontSize', 6);
                set(xl, 'FontSize', xlFontSize);
                set(yl, 'FontSize', xlFontSize);

                text(110,yt,[' f = ' num2str(fi,'%.2f') ' Hz'],'edgecolor','k','backgroundcolor','w')
                jeffsavefig_r300(fat2,[SaveDirString '/Region' num2str(Options.JeffRegionFlag) '/Figures/' num2str(Evi,'%d') '_Atten2' RegString '.png']);
                if Options.CloseFigs
                    close(fat2)
                end
            end
            
            
        end
        
    end
    hai=[2 3 5 6 8 9];
    set(ha(1:6),'XTickLabel',''); %set(ha(hai),'YTickLabel','')
    set(ha,'XMinorGrid','off')
    set(ha,'XMinorTick','off')
%     axes(ha(1));
%     text(75,1.8e-2,'(b)','fontsize',20);
%         yAX = get(gca,'YAxis');
%         set(yAX,'FontSize', 6);
    jeffsavefig_r300(fhdist,[SaveDirString '/Region' num2str(Options.JeffRegionFlag) '/Figures/' num2str(Evi,'%d') '_Atten' RegString '.png']);
    if Options.CloseFigs
        close(fhdist)
    end

    %% now, for this event, fit the freq dependence of Q
    Qposi=Qall(ii,:)>0 & QallPlus(ii,:)>0  & QallMinus(ii,:)>0 & QallPlus(ii,:)<10000; % indices of positive Q, don't use negative instances, remove ridiculous ones due to c being very close to 0   
    if SpecialEventFlag==7 || SpecialEventFlag==10 
        Qposi=Qall(ii,:)>0 & QallMinus(ii,:)>0 & QallPlus(ii,:)<10000;
    end
    
    if sum(Qposi)<3
        disp(['Skipping Event ' num2str(Evi) ', ' EvNamei ', Bad Q Calc'])
        continue
    elseif sum(Qposi)<size(FASi,2)
        qnadd=logical(ones(size(Qposi))-Qposi);
        Qall(ii,qnadd)=nan;
        QallPlus(ii,qnadd)=nan;
        QallMinus(ii,qnadd)=nan;
    end

    Y2=log(Qall(ii,Qposi)).';
    beta2=[100;0.5]; % initial values for the fit
    [beta3,beta3res,jaco3,covb3,mse3]=nlinfit(FREQS(Qposi).',Y2,@QFIT,beta2,options); % linear fit
%     [beta4,beta4res,jaco4,covb4,mse4]=nlinfit(FREQS(Qposi).',Y2,@QFIT2,beta2,options); % exponential fit?
    xdomax=1.0843*max(FREQS(Qposi));
    xdo=logspace(log10(1.0),log10(xdomax),100); % range of freq values
    [ydo2,ydelta2] = nlpredci(@QFIT,xdo.',beta3,beta3res,'covar',covb3); % prediction for these freq, with 95% CI
%     [ydo3,ydelta3] = nlpredci(@QFIT2,xdo.',beta4,beta4res,'covar',covb4); % prediction for these freq, with 95% CI, exponential fit

    % save Qo and eta
    Qo(ii,1)=beta3(1); 
    Qo(ii,2)=sqrt(covb3(1,1)); % the sigma
    eta(ii,1)=beta3(2);
    eta(ii,2)=sqrt(covb3(2,2)); % the sigma

    if Options.OtherEASAdj % do the fit for EAS_raw and EAS_Site and EAS_rad
        
        % for site adjusted data, EAS_site
        QposiS=QallS(ii,:)>0 & QallPlusS(ii,:)>0  & QallMinusS(ii,:)>0 & QallPlusS(ii,:)<10000; 
        if SpecialEventFlag==7 || SpecialEventFlag==10 
            QposiS=QallS(ii,:)>0 & QallMinusS(ii,:)>0 & QallPlusS(ii,:)<10000;
        end
        if sum(QposiS)<3
%             disp(['Skipping Event ' num2str(Evi) ', ' EvNamei ', Bad Q Calc'])
            continue
        elseif sum(QposiS)<size(FASi,2)
            qnadd=logical(ones(size(QposiS))-QposiS);
            QallS(ii,qnadd)=nan;
            QallPlusS(ii,qnadd)=nan;
            QallMinusS(ii,qnadd)=nan;
        end

        YS=log(QallS(ii,QposiS)).';
        beta2=[100;0.5]; % initial values for the fit
        [beta3S,beta3resS,jaco3S,covb3S,mse3S]=nlinfit(FREQS(QposiS).',YS,@QFIT,beta2,options); % linear fit
        % save Qo and eta
        QoS(ii,1)=beta3S(1); 
        QoS(ii,2)=sqrt(covb3S(1,1)); % the sigma
        etaS(ii,1)=beta3S(2);
        etaS(ii,2)=sqrt(covb3S(2,2)); % the sigma
        
        % for EAS_rad
        QposiR=QallR(ii,:)>0 & QallPlusR(ii,:)>0  & QallMinusR(ii,:)>0 & QallPlusR(ii,:)<10000; 
        if SpecialEventFlag==7 || SpecialEventFlag==10 
            QposiR=QallR(ii,:)>0 & QallMinusR(ii,:)>0 & QallPlusR(ii,:)<10000;
        end
        if sum(QposiR)<3
%             disp(['Skipping Event ' num2str(Evi) ', ' EvNamei ', Bad Q Calc'])
            continue
        elseif sum(QposiR)<size(FASi,2)
            qnadd=logical(ones(size(QposiR))-QposiR);
            QallR(ii,qnadd)=nan;
            QallPlusR(ii,qnadd)=nan;
            QallMinusR(ii,qnadd)=nan;
        end

        YR=log(QallR(ii,QposiR)).';
        beta2=[100;0.5]; % initial values for the fit
        [beta3R,beta3resR,jaco3R,covb3R,mse3R]=nlinfit(FREQS(QposiR).',YR,@QFIT,beta2,options); % linear fit
        % save Qo and eta
        QoR(ii,1)=beta3R(1); 
        QoR(ii,2)=sqrt(covb3R(1,1)); % the sigma
        etaR(ii,1)=beta3R(2);
        etaR(ii,2)=sqrt(covb3R(2,2)); % the sigma
        
        % for EAS_raw
        QposiRaw=QallRaw(ii,:)>0 & QallPlusRaw(ii,:)>0  & QallMinusRaw(ii,:)>0 & QallPlusRaw(ii,:)<10000; 
        if SpecialEventFlag==7 || SpecialEventFlag==10 
            QposiRaw=QallRaw(ii,:)>0 & QallMinusRaw(ii,:)>0 & QallPlusRaw(ii,:)<10000;
        end
        if sum(QposiRaw)<3
%             disp(['Skipping Event ' num2str(Evi) ', ' EvNamei ', Bad Q Calc'])
            continue
        elseif sum(QposiRaw)<size(FASi,2)
            qnadd=logical(ones(size(QposiRaw))-QposiRaw);
            QallRaw(ii,qnadd)=nan;
            QallPlusRaw(ii,qnadd)=nan;
            QallMinusRaw(ii,qnadd)=nan;
        end

        YRaw=log(QallRaw(ii,QposiRaw)).';
        beta2=[100;0.5]; % initial values for the fit
        [beta3Raw,beta3resRaw,jaco3Raw,covb3Raw,mse3Raw]=nlinfit(FREQS(QposiRaw).',YRaw,@QFIT,beta2,options); % linear fit
        % save Qo and eta
        QoRaw(ii,1)=beta3Raw(1); 
        QoRaw(ii,2)=sqrt(covb3Raw(1,1)); % the sigma
        etaRaw(ii,1)=beta3Raw(2);
        etaRaw(ii,2)=sqrt(covb3Raw(2,2)); % the sigma
                
    end

    fhQ=figure(1000+ii); set(fhQ,'position',[12   563   759   360]);%[12   580   555   343]);
    loglog(FREQS(Qposi),Qall(ii,Qposi),'ko','markerfacecolor',rgb('indianred'),'markersize',8); hold on
    loglog(FREQS(Qposi),QallPlus(ii,Qposi),'kv','markerfacecolor',rgb('indianred'),'markersize',8);
    loglog(FREQS(Qposi),QallMinus(ii,Qposi),'k^','markerfacecolor',rgb('indianred'),'markersize',8);
    loglog(FREQS(Qposi),Qall(ii,Qposi),'ko','markerfacecolor',rgb('indianred'),'markersize',8); % replot for order
    grid on; hold on
    loglog(xdo,exp(ydo2),'-','color',rgb('black'),'linewidth',1)
    loglog(xdo,exp(ydo2+ydelta2),'--','color',rgb('black'),'linewidth',1)
    loglog(xdo,exp(ydo2-ydelta2),'--','color',rgb('black'),'linewidth',1)
    xlabel('\itf\rm (Hz)')
    ylabel(['Apparent \itQ(f)'])
    axis([.8 20 100 4e4])
    title(['M' num2str(Mi,'%.1f') ' ' strrep(EvNamei,'_',' ') ', EQID ' num2str(Evi,'%d') ', Region: ' RegionString ])
    tstr=char(['Q_0 = ' num2str(round(Qo(ii,1)),'%.d') ' \pm ' num2str(round(Qo(ii,2)),'%.d')],...
        ['\eta = ' num2str(eta(ii,1),'%.2f') ' \pm ' num2str(eta(ii,2),'%.2f')]);

    text(1,17000,tstr,'edgecolor','k','backgroundcolor','w','fontsize',14)
    text(.67,5e4,'(d)','fontsize',20);
    
    removewhitespace(gca,1);
    jeffsavefig_r300(fhQ,[SaveDirString '/Region' num2str(Options.JeffRegionFlag) '/Figures/' num2str(Evi,'%d') '_Q' RegString '.png']);
    if Options.CloseFigs
        close(fhQ)
    end

    % plot the frequency dependence of the anelastic attenutation coeff
    fhC=figure(10000+ii); set(fhC,'position',[12   563   759   360]);%[12   580   555   343]);
    semilogx(FREQS(Qposi),Call(ii,Qposi),'ko','markerfacecolor',rgb('forestgreen'),'markersize',8); hold on; grid on
    loglog(FREQS(Qposi),CallPlus(ii,Qposi),'kv','markerfacecolor',rgb('forestgreen'),'markersize',8);
    loglog(FREQS(Qposi),CallMinus(ii,Qposi),'k^','markerfacecolor',rgb('forestgreen'),'markersize',8);
    xlabel('\itf\rm (Hz)')
    ylabel(['Anelastic Attenuation Coefficient, \itc'])
    axis([.8 20 -.02 0])
    title(['M' num2str(Mi,'%.1f') ' ' strrep(EvNamei,'_',' ') ', EQID ' num2str(Evi,'%d') ', Region: ' RegionString ]) 
    text(.63,.0009,'(c)','fontsize',20);
    
    removewhitespace(gca,1);
    jeffsavefig_r300(fhC,[SaveDirString '/Region' num2str(Options.JeffRegionFlag) '/Figures/' num2str(Evi,'%d') '_C' RegString '.png']);
    if Options.CloseFigs
        close(fhC)
    end
        
    % save some data for this info in a mat file
    
    Ev.Qo=Qo(ii,:);
    Ev.eta=eta(ii,:);
    Ev.f=FREQS;
    Ev.C=Call(ii,:);
    Ev.CSig=.5*(CallPlus(ii,:)-CallMinus(ii,:));
    Ev.Q=Qall(ii,:);
    Ev.QP=QallPlus(ii,:);
    Ev.QM=QallMinus(ii,:);
    Ev.M=Mi;
    Ev.Ritrack=Ritrack;
    Ev.Evname=EvNamei;
    Ev.Evi=Evi;
    Ev.Rmax=Rmax;
    Ev.Rmin=Rmin;
    Ev.RadType=Options.RadType;
    Ev.JRegion=Options.JeffRegionFlag;
    
    
    if Options.OtherEASAdj % save results for EAS_raw and EAS_Site
        Ev.QoS=QoS(ii,:);
        Ev.etaS=etaS(ii,:);
        Ev.QS=QallS(ii,:);
        Ev.CSigS=.5*(CallPlusS(ii,:)-CallMinusS(ii,:));
        
        Ev.QoR=QoR(ii,:);
        Ev.etaR=etaR(ii,:);
        Ev.QR=QallR(ii,:);
        Ev.CSigR=.5*(CallPlusR(ii,:)-CallMinusR(ii,:));
        
        Ev.QoRaw=QoRaw(ii,:);
        Ev.etaRaw=etaRaw(ii,:);
        Ev.QRaw=QallRaw(ii,:);
        Ev.CSigRaw=.5*(CallPlusRaw(ii,:)-CallMinusRaw(ii,:));
        
        Ev.ysig=YSIG(ii,:);
        Ev.ysigR=YSIGR(ii,:);
        Ev.ysigS=YSIGS(ii,:);
        Ev.ysigRaw=YSIGRAW(ii,:);
    end
    
    save([SaveDirString '/Region' num2str(Options.JeffRegionFlag) '/' num2str(Evi,'%d') '_Results' RegString '.mat'],'Ev')
 
    
    %% finish the map
    % Cramer (2014) boundaries
    figure(1)
    pcol='amethyst';
    pp1=m_plot(PEERceus(:,1),PEERceus(:,2),'color',rgb(pcol),'linewidth',2);
    pp2=m_plot(PEERgulf(:,1),PEERgulf(:,2),'color',rgb(pcol),'linewidth',2);
    pp3=m_plot(PEERapp(:,1),PEERapp(:,2),'color',rgb(pcol),'linewidth',2);
    pp4=m_plot(PEERacp(:,1),PEERacp(:,2),'color',rgb(pcol),'linewidth',2);

    lh=legend([pp1,prad],char('NGA-East Regions Boundaries','     (Dreiling et al, 2014)'),...
        'S Radiation Pattern, f=1.5 Hz');
    set(lh,'fontsize',10)
    % grid 
    m_grid('box','fancy','tickdir','in');
    m_ruler([-.075 .1],.001,'tickdir','out','ticklen',[.007 .007],'fontsize',10);
    title(char(['M' num2str(Mi,'%.1f') ' ' strrep(EvNamei,'_',' ') ', EQID ' num2str(Evi,'%d') ', Region: ' RegionString  ],...
        ['                 [\phi=' num2str(Strikei) ' \delta=' num2str(Dipi) ' \lambda=' num2str(Rakei)  ']']))
    ax2=gca;
    ax2.OuterPosition=[-.055 -.08 1.13 1.13];
    text(-.45,.35,'(a)','fontsize',20);

    jeffsavefig_r300(f1,[SaveDirString '/Region' num2str(Options.JeffRegionFlag) '/Figures/' num2str(Evi,'%d') '_Map' RegString '.png']);
    if Options.CloseFigs
        close(f1)
    end
 
end

disp('finished with script')

