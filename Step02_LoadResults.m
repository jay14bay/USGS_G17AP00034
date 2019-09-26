%% script for USGS Grant G17AP00034
% Jeff Bayless, AECOM, jeff.bayless@aecom.com
% 26 September 2019

% this aggregates the results from script "Step01_MainScriptForFAS.m"

close all; clearvars -except PNUM2 PTXT2 PRAW2; clc;

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
% select region
Options.JeffRegionFlag = 3 ; % this is a flag for the region used in selecting the data to analyze.
        % 1 = Gulf Coast, source/site crosses the gulf coast boundary
        % 2 = cena
        % 3 = appalachian province
        % 4 = Atlantic Coastal Plain
        % 5 = all data, irrespective of region
        % 11 = gulf coast; strictly path 1 only

%% Event selection
if Options.JeffRegionFlag==11 % Gulf Coast only (strict)
    GCEQID=   [46	47	48	49	66	67	73	74	76	80	81	90	91	92];
    UseForAvg=[1     1   1   1   1   1   1   1   1   1   1   1   1   1];
    EvUnique=GCEQID;
    tsr='(a)';
    RegionString='Region 1: Gulf Coast';
elseif Options.JeffRegionFlag==2 % CENA
    CENAID=   [16 25 29 30 33 44 46 47 55 56 58 60 61 66 67 73 74 75 76 80 81 83 90 91]; % didnt use: 21 32 34 35 37 51 85 57
    UseForAvg=[1   1  1  1  1  1  0  1  1  1  1  0  0  1  1  1  1  1  1  1  1  1  1  1];
    EvUnique=CENAID;
    tsr='(b)';
    RegionString='Region 2: CNA';
elseif Options.JeffRegionFlag==3 % Appalachian
    APPID=    [16 26 51 68 82 85 88 89 ];
    UseForAvg=[ 1  1  1  1  1  1  1  1 ];
    EvUnique=APPID;
    tsr='(c)';
    RegionString='Region 3: Appalachian Province';
elseif Options.JeffRegionFlag==4 % atlantic coast
    EvUnique=[21 88];
    UseForAvg=[1 1];
elseif Options.JeffRegionFlag==5 % all data

end
RegString=['Region' num2str(Options.JeffRegionFlag) ];
lstr=[]; cc=0;

SigDif_Raw_RadSite=nan(sum(UseForAvg),10);
SigDif_Raw_Site=nan(sum(UseForAvg),10);
SigDif_Raw_Rad=nan(sum(UseForAvg),10);
CSE_RadSite=nan(sum(UseForAvg),10);
CSE_Site=nan(sum(UseForAvg),10);
CSE_Rad=nan(sum(UseForAvg),10);

for ii=1:length(EvUnique) 
    Evi=EvUnique(ii);
    load([SaveDirString '/' RegString '/' num2str(Evi,'%d') '_Results_' RegString '.mat'])

    figure(1)
    nr=length(Ev.Ritrack);
    semilogx(Ev.Ritrack,repmat(Ev.M,nr),'ko','markerfacecolor','b'); hold on
    xlim([100 600])
    
    figure(2)
        subplot(300,1,1:80)
        semilogy(Ev.M,Ev.Qo(1),'ko','markerfacecolor','k'); hold on
        grid on
        ylabel('\itQ_0')
        ylim([1e2 1e4])
        xlabel('\bfM')
    subplot(300,1,101:180)
        plot(Ev.M,Ev.eta(1),'ko','markerfacecolor','k'); hold on
        grid on
        xlabel('\bfM')
        ylabel('\it\eta')
        ylim([-0.5 1.2])
    subplot(300,1,201:280)
        plot(Ev.M,Ev.C(6),'ko','markerfacecolor','k'); hold on % 5 hz
        grid on
        xlabel('\bfM')
        ylabel('\itc')
        ylim([-0.01 0])
    
    SaveTable(ii,:)=[Evi Ev.Qo Ev.eta];
    Qtable(ii,:)=Ev.Q;
    QPtable(ii,:)=Ev.QP;
    QMtable(ii,:)=Ev.QM;
    
    if UseForAvg(ii)==1
        cc=cc+1;
        lstr{cc}=['EQID ' num2str(Evi)];
    
        % estimates from EASraw and EASsite
        SaveTableS(cc,:)=[Evi Ev.QoS Ev.etaS];
        SaveTableR(cc,:)=[Evi Ev.QoR Ev.etaR];
        SaveTableRaw(cc,:)=[Evi Ev.QoRaw Ev.etaRaw];
    
        SigDif_Raw_RadSite(cc,1:length(Ev.ysigRaw))=[Ev.ysigRaw-Ev.ysig];
        SigDif_Raw_Site(cc,1:length(Ev.ysigRaw))=[Ev.ysigRaw-Ev.ysigS];
        SigDif_Raw_Rad(cc,1:length(Ev.ysigRaw))=[Ev.ysigRaw-Ev.ysigR];

        CSE_RadSite(cc,1:length(Ev.ysigRaw))=[Ev.CSigRaw-Ev.CSig];
        CSE_Site(cc,1:length(Ev.ysigRaw))=[Ev.CSigRaw-Ev.CSigS];
        CSE_Rad(cc,1:length(Ev.ysigRaw))=[Ev.CSigRaw-Ev.CSigR];
    
    end
    
end
f=Ev.f;

SIGD1=mean(SigDif_Raw_RadSite,1);
SIGD2=mean(SigDif_Raw_Site,1);
SIGD3=mean(SigDif_Raw_Rad,1);

figure; set(gcf,'position',[553    66   560   919]);
subplot(28,1,1:8)
    semilogx(f,SigDif_Raw_RadSite,'o-'); hold on
    semilogx(f,SIGD1,'k-','linewidth',2)
    semilogx(f,zeros(size(f)),'k-')
    title('\sigma Reduction: \itEAS_{raw} - EAS_{RadSite}')
    ylim([-.25 .25])
    grid on
    xlabel('frequency (Hz)')
subplot(28,1,11:18)
    semilogx(f,SigDif_Raw_Site,'o-'); hold on
    semilogx(f,SIGD2,'k-','linewidth',2)
    semilogx(f,zeros(size(f)),'k-')
    title('\sigma Reduction: \itEAS_{raw} - EAS_{Site}')
    ylim([-.25 .25])
    grid on
    xlabel('frequency (Hz)')
subplot(28,1,21:28)
    semilogx(f,SigDif_Raw_Rad,'o-'); hold on
    semilogx(f,SIGD3,'k-','linewidth',2)
    semilogx(f,zeros(size(f)),'k-')
    title('\sigma Reduction: \itEAS_{raw} - EAS_{Rad}')
    ylim([-.25 .25])
    grid on
    xlabel('frequency (Hz)')

CSE1=mean(CSE_RadSite,1);
CSE2=mean(CSE_Site,1);
CSE3=mean(CSE_Rad,1);
figure; set(gcf,'position',[553    66   560   919]);
subplot(28,1,1:8)
    semilogx(f,CSE_RadSite,'o-'); hold on
    semilogx(f,CSE1,'k-','linewidth',2)
    semilogx(f,zeros(size(f)),'k-')
    title('se_c Reduction: \itEAS_{raw} - EAS_{RadSite}')
    ylim([-.5e-3 .5e-3])
    grid on
    xlabel('frequency (Hz)')
subplot(28,1,11:18)
    semilogx(f,CSE_Site,'o-'); hold on
    semilogx(f,CSE2,'k-','linewidth',2)
    semilogx(f,zeros(size(f)),'k-')
    title('se_c Reduction: \itEAS_{raw} - EAS_{Site}')
    ylim([-.5e-3 .5e-3])
    grid on
    xlabel('frequency (Hz)')
subplot(28,1,21:28)
    semilogx(f,CSE_Rad,'o-'); hold on
    semilogx(f,CSE3,'k-','linewidth',2)
    semilogx(f,zeros(size(f)),'k-')
    title('se_c Reduction: \itEAS_{raw} - EAS_{Rad}')
    ylim([-.5e-3 .5e-3])
    grid on
    xlabel('frequency (Hz)')

lstr{cc+1}='Mean';
figure; set(gcf,'position',[541   202   598   708]);
subplot(17,1,1:8)
    semilogx(f,SigDif_Raw_RadSite,'-'); hold on
    semilogx(f,SIGD1,'k-','linewidth',2)
    semilogx(f,zeros(size(f)),'k-')
%     title('\sigma Reduction: \itEAS_{raw} - EAS_{RadSite}')
    title('\sigma Reduction')
    ylabel('ln units')
    ylim([-.2 .2].*2)
    grid on
    lh=legend(lstr,'location','eastoutside')
%     xlabel('frequency (Hz)')
subplot(17,1,10:17)
    semilogx(f,CSE_RadSite,'-'); hold on
    semilogx(f,CSE1,'k-','linewidth',2)
    semilogx(f,zeros(size(f)),'k-')
%     title('se_c Reduction: \itEAS_{raw} - EAS_{RadSite}')
    title('se_c Reduction')
    ylim([-.6e-3 .6e-3].*2)
    grid on
    xlabel('frequency (Hz)')
    lh=legend(lstr,'location','eastoutside')
    
jeffsavefig_r300(gcf,[SaveDirString '/Region' num2str(Options.JeffRegionFlag) '/SIGMA_' RegString '.png']);

%% Re-calculate Q(f) from event-averaged Q

% fit the mean Q(f) again (unweighted)
Q0=Qtable(logical(UseForAvg),:); % remove the ones we want to skip
Qneg=Q0<=0; %negative entries
Q0(Qneg)=nan; %makes the negative entries nan
Q1=nanmean(log(Q0),1); % average of events, excluding nan's
Y2=Q1.'; Qall=exp(Y2);

QP0=QPtable(logical(UseForAvg),:); % remove the ones we want to skip
Qneg=QP0<=0; %negative entries
QP0(Qneg)=nan; %makes the negative entries nan
QP1=nanmean(log(QP0),1); % average of events
QPall=exp(QP1.');

QM0=QMtable(logical(UseForAvg),:); % remove the ones we want to skip
Qneg=QM0<=0; %negative entries
QM0(Qneg)=nan; %makes the negative entries nan
QM1=nanmean(log(QM0),1); % average of events
QMall=exp(QM1.');

beta2=[500;0.5]; % initial values for the fit
options = statset('TolX',1e-15,'TolFun',1e-15,'robust','on');
options.RobustWgtFun = 'huber'; %'bisquare';
[beta3,beta3res,jaco3,covb3,mse3]=nlinfit(f.',Y2,@QFIT,beta2,options); % nonlinear fit using Cramer (2014) eqn 1
%     [beta4,beta4res,jaco4,covb4,mse4]=nlinfit(f.',Y2,@QFIT2,beta2,options); % eponential fit 
xdomax=1.0843*max(f);
xdo=logspace(log10(1.0),log10(xdomax),100); % range of freq values
[ydo2,ydelta2] = nlpredci(@QFIT,xdo.',beta3,beta3res,'covar',covb3); % prediction for these freq, with 95% CI
%     [ydo3,ydelta3] = nlpredci(@QFIT2,xdo.',beta4,beta4res,'covar',covb4); % prediction for these freq, with 95% CI, exponential fit

% save Qo and eta
Qo(1,1)=beta3(1); 
Qo(1,2)=sqrt(covb3(1,1)); % the sigma
eta(1,1)=beta3(2);
eta(1,2)=sqrt(covb3(2,2)); % the sigma

fhQ=figure(1000+ii); set(fhQ,'position',[12   563   901   360]);%[12   580   555   343]);
loglog(f,Qall,'ko','markerfacecolor',rgb('indianred'),'markersize',8); hold on
loglog(f,QPall,'kv','markerfacecolor',rgb('indianred'),'markersize',8);
loglog(f,QMall,'k^','markerfacecolor',rgb('indianred'),'markersize',8);
loglog(f,Qall,'ko','markerfacecolor',rgb('indianred'),'markersize',8); % replot for order
grid on; hold on
loglog(xdo,exp(ydo2),'-','color',rgb('black'),'linewidth',1)
loglog(xdo,exp(ydo2+ydelta2),'--','color',rgb('black'),'linewidth',1)
loglog(xdo,exp(ydo2-ydelta2),'--','color',rgb('black'),'linewidth',1)
xlabel('\itf\rm (Hz)')
ylabel(['Apparent \itQ(f)'])
axis([.8 20 100 4e4])
% title(['M' num2str(Mi,'%.1f') ' ' strrep(EvNamei,'_',' ') ', EQID ' num2str(Evi,'%d') ', Region: ' RegionString ])
tstr=char(['Q_0 = ' num2str(round(Qo(1,1)),'%.d') ' \pm ' num2str(round(Qo(1,2)),'%.d')],...
    ['\eta = ' num2str(eta(1,1),'%.2f') ' \pm ' num2str(eta(1,2),'%.2f')]);
%     tstr2=char(['a = ' num2str(round(aexp(ii,1)),'%.d') ' \pm ' num2str(round(aexp(ii,2)),'%.d')],...
%         ['b = ' num2str(bexp(ii,1),'%.2f') ' \pm ' num2str(bexp(ii,2),'%.1e')]);
text(1,17000,tstr,'edgecolor','k','backgroundcolor','w','fontsize',14)
text(.7,5e4,tsr,'fontsize',20);
title(RegionString)

removewhitespace(gca,1);
jeffsavefig_r300(fhQ,[SaveDirString '/Region' num2str(Options.JeffRegionFlag) '/AverageQ_' RegString '.png']);

disp('finished')
    