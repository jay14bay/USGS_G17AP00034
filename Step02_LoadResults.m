close all; clearvars -except PNUM2 PTXT2 PRAW2; clc;

%% select region
Options.JeffRegionFlag = 11 ; % this is a flag for the region used in selecting the data to analyze.
        % 1 = Gulf Coast, source/site crosses the gulf coast boundary
        % 2 = cena
        % 3 = appalachian province
        % 4 = Atlantic Coastal Plain
        % 5 = all data, irrespective of region
        % 11 = gulf coast; strictly path 1 only

%% Event selection
if Options.JeffRegionFlag==11 % Gulf Coast only (strict)
    GCEQID=   [46	47	48	49	58	66	67	73	74	76	80	81	90	91	92];
    UseForAvg=[1     1   1   1   1   1   1   1   1   1   1   1   1   1   0];
    EvUnique=GCEQID;
elseif Options.JeffRegionFlag==2 % CENA
    CENAID=   [16 25 29 30 33 34 35 44 46 47 48 55 56 57 58 60 61 66 67 73 74 75 76 80 81 83 85 90 91]; 
    UseForAvg=[1   1  1  1  1  1  0  1  0  1  0  1  1  1  1  0  1  1  0  1  1  1  0  0  0  1  0  1  0];
    EvUnique=CENAID;
elseif Options.JeffRegionFlag==3 % Appalachian
    
elseif Options.JeffRegionFlag==4 % atlantic coast
    
elseif Options.JeffRegionFlag==5 % all data

end
RegString=['Region' num2str(Options.JeffRegionFlag)];

for ii=1:length(EvUnique) 
    Evi=EvUnique(ii);
    load([RegString '/' num2str(Evi,'%d') '_Results_' RegString '.mat'])
    
    figure(1)
    nr=length(Ev.Ritrack);
    semilogx(Ev.Ritrack,repmat(Ev.M,nr),'ko','markerfacecolor','b'); hold on
    xlim([100 600])
    
    figure(2)
    semilogy(Ev.M,Ev.Qo,'ko','markerfacecolor','k'); hold on
    
    MassiveTable(ii,:)=[Evi Ev.Qo Ev.eta];
    Qtable(ii,:)=Ev.Q;
    QPtable(ii,:)=Ev.QP;
    QMtable(ii,:)=Ev.QM;
    
    
end


%% Re-calculate Q(f) from event-averaged Q

f=Ev.f;

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

fhQ=figure(1000+ii); set(fhQ,'position',[12   563   759   360]);%[12   580   555   343]);
loglog(f,Qall,'ko','markerfacecolor',rgb('indianred'),'markersize',8); hold on
loglog(f,QPall,'kv','markerfacecolor',rgb('indianred'),'markersize',8);
loglog(f,QMall,'k^','markerfacecolor',rgb('indianred'),'markersize',8);
loglog(f,Qall,'ko','markerfacecolor',rgb('indianred'),'markersize',8); % replot for order
grid on; hold on
loglog(xdo,exp(ydo2),'-','color',rgb('black'),'linewidth',1)
loglog(xdo,exp(ydo2+ydelta2),'--','color',rgb('black'),'linewidth',1)
loglog(xdo,exp(ydo2-ydelta2),'--','color',rgb('black'),'linewidth',1)
%         loglog(xdo,exp(ydo3),'-','color',rgb('black'),'linewidth',1.5)
%         loglog(xdo,exp(ydo3+ydelta3),'--','color',rgb('black'),'linewidth',1.5)
%         loglog(xdo,exp(ydo3-ydelta3),'--','color',rgb('black'),'linewidth',1.5)
xlabel('\itf\rm (Hz)')
ylabel(['Apparent \itQ(f)'])
axis([.8 20 100 4e4])
% title(['M' num2str(Mi,'%.1f') ' ' strrep(EvNamei,'_',' ') ', EQID ' num2str(Evi,'%d') ', Region: ' RegionString ])
tstr=char(['Q_0 = ' num2str(round(Qo(1,1)),'%.d') ' \pm ' num2str(round(Qo(1,2)),'%.d')],...
    ['\eta = ' num2str(eta(1,1),'%.2f') ' \pm ' num2str(eta(1,2),'%.2f')]);
%     tstr2=char(['a = ' num2str(round(aexp(ii,1)),'%.d') ' \pm ' num2str(round(aexp(ii,2)),'%.d')],...
%         ['b = ' num2str(bexp(ii,1),'%.2f') ' \pm ' num2str(bexp(ii,2),'%.1e')]);
text(1,17000,tstr,'edgecolor','k','backgroundcolor','w','fontsize',14)
% text(.67,5e4,'(d)','fontsize',20);

removewhitespace(gca,1);
jeffsavefig_r300(fhQ,['Region' num2str(Options.JeffRegionFlag) '/AverageQ_' RegString '.png']);

%% compare with other models




disp('finished')
    