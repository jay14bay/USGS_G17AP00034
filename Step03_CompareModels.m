%% script for USGS Grant G17AP00034
% Jeff Bayless, AECOM, jeff.bayless@aecom.com
% 26 September 2019

% this plots the results from script "Step02_LoadResults.m"

close all; clearvars -except PNUM2 PTXT2 PRAW2; clc;

f=[1.0000    1.5100    2.1400    2.8100    3.8000    5.0100    6.7600    9.1200   12.5900   16.6000];
xdomax=1.0843*max(f);
xdo=logspace(log10(1.0),log10(xdomax),100); % range of freq values
yma=5e3;

posnarrow=[12   563   901   360];
postall=[2088          35         937         853]

%% for Gulf Coast first

%mine
Q=278;
eta=0.6;

fh1=figure(1); set(fh1,'position',posnarrow);
loglog(xdo,Q.*xdo.^eta,'--','color',rgb('black'),'linewidth',2.5); grid on; hold on
lstr{1}='Gulf Coast; This Study [b=-0.5]';

%Cramer12
loglog(xdo,270.*xdo.^.75,'--','color',rgb('dodgerblue'),'linewidth',2);
lstr=[lstr 'Gulf Coast; Cramer (2012) [b=-0.5]'];

%Silva 02
loglog(xdo,351.*xdo.^.84,'--','color',rgb('indianred'),'linewidth',2);
lstr=[lstr 'Gulf Coast; Silva et al. (2002) [b=-0.55]' ];


xlabel('\itf\rm (Hz)')
ylabel(['\itQ(f)'])
axis([.8 20 100 yma])
% title('Region 1: Gulf Coast')
text(.67,4.5e3,'(a)','fontsize',20);
 
lh=legend(lstr,'location','northwest'); set(lh,'fontsize',13)
ah=gca; set(ah,'fontsize',13)

removewhitespace(gca,1);
jeffsavefig_r300(fh1,['Q_Comparison_R1.png']);


%% for CNA
clear lstr
%mine
Q=465;
eta=0.56;

fh1=figure(10); set(fh1,'position',posnarrow);%[12   580   555   343]);
loglog(xdo,Q.*xdo.^eta,'-','color',rgb('black'),'linewidth',2.5); grid on; hold on
lstr{1}='CNA; This Study [b=-0.5]';

%Cramer 2012
loglog(xdo,550.*xdo.^.69,'-','color',rgb('dodgerblue'),'linewidth',2);
lstr=[lstr 'CNA; Cramer (2012) [b=-0.5]'];

%Erickson 2004
loglog(xdo,640.*xdo.^.34,'-','color',rgb('indianred'),'linewidth',2);
lstr=[lstr 'CNA; Erickson et al. (2004) [b=-0.5]'];


%% for ENA

%Cramer 2012
loglog(xdo,550.*xdo.^.6,'-.','color',rgb('purple'),'linewidth',2);
lstr=[lstr 'ENA; Cramer (2012) [b=-0.5]'];

%AB95
loglog(xdo,680.*xdo.^.36,'-.','color',rgb('forestgreen'),'linewidth',2);
lstr=[lstr 'ENA; Atkinson and Boore (1995) [b=-0.5]'];

%AB14
loglog(xdo,525.*xdo.^.45,'-.','color',rgb('cyan'),'linewidth',2);
lstr=[lstr 'ENA; Atkinson and Boore (2014) [b=-0.5]'];

%BS11
loglog(xdo,410.*xdo.^.5,'-.','color',rgb('magenta'),'linewidth',2);
lstr=[lstr 'ENA; Boatwright and Seekins (2011) [b=-0.5]'];

%Atkinson 2004
loglog(xdo,max(repmat(1000,size(xdo)),893.*xdo.^.32),'-.','color',rgb('orange'),'linewidth',2);
lstr=[lstr 'ENA; Atkinson (2004) [b=-0.5]'];

%Erickson 2004
loglog(xdo,650.*xdo.^.36,'-.','color',rgb('brown'),'linewidth',2);
lstr=[lstr 'ENA; Erickson et al. (2004) [b=-0.5]'];

%replot mine on top
loglog(xdo,Q.*xdo.^eta,'-','color',rgb('black'),'linewidth',2.5);


xlabel('\itf\rm (Hz)')
ylabel(['\itQ(f)'])
axis([.8 20 100 yma])
% title('Region 1: Gulf Coast')
text(.67,4.5e3,'(b)','fontsize',20);
 
lh=legend(lstr,'location','southeast'); set(lh,'fontsize',13)
ah=gca; set(ah,'fontsize',13)

removewhitespace(gca,1);
jeffsavefig_r300(fh1,['Q_Comparison_R2.png']);

%% for Appalachian
clear lstr
%mine
Q=451;
eta=0.55;

fh1=figure(100); set(fh1,'position',posnarrow);%[12   580   555   343]);
loglog(xdo,Q.*xdo.^eta,':','color',rgb('black'),'linewidth',2.5); grid on; hold on
lstr{1}='Appalachian; This Study [b=-0.5]';

%Shi 96
loglog(xdo,573.5.*xdo.^.465,':','color',rgb('dodgerblue'),'linewidth',2);
lstr=[lstr 'Appalachian; Shi et al. (1996) [b=n/a]'];

xlabel('\itf\rm (Hz)')
ylabel(['\itQ(f)'])
axis([.8 20 100 yma])
% title('Region 1: Gulf Coast')
text(.67,4.5e3,'(c)','fontsize',20);
 
lh=legend(lstr,'location','northwest'); set(lh,'fontsize',13)
ah=gca; set(ah,'fontsize',13)

removewhitespace(gca,1);
jeffsavefig_r300(fh1,['Q_Comparison_R3.png']);


% 
% 
% 
% %% all on one
% 
% %% for Gulf Coast first
% 
% %mine
% Q=304;
% eta=0.6;
% 
% fh1=figure(1); set(fh1,'position',[2088          35         937         853]);%[12   580   555   343]);
% loglog(xdo,Q.*xdo.^eta,'--','color',rgb('black'),'linewidth',2.5); grid on; hold on
% lstr{1}='Gulf Coast; This Study [b=-0.5]';
% 
% %Cramer12
% loglog(xdo,270.*xdo.^.75,'--','color',rgb('dodgerblue'),'linewidth',2);
% lstr=[lstr 'Gulf Coast; Cramer (2012) [b=-0.5]'];
% 
% %Silva 02
% loglog(xdo,351.*xdo.^.84,'--','color',rgb('indianred'),'linewidth',2);
% lstr=[lstr 'Gulf Coast; Silva et al. (2002) [b=-0.55]' ];
% 
% plot(99,9999,'wo','markerfacecolor','w'); lstr=[lstr newline];
% 
% %% for CNA
% 
% %mine
% Q=433;
% eta=0.58;
% loglog(xdo,Q.*xdo.^eta,'-','color',rgb('black'),'linewidth',2.5);
% lstr=[lstr 'CNA; This Study [b=-0.5]'];
% 
% %Cramer 2012
% loglog(xdo,550.*xdo.^.69,'-','color',rgb('dodgerblue'),'linewidth',2);
% lstr=[lstr 'CNA; Cramer (2012) [b=-0.5]'];
% 
% %Erickson 2004
% loglog(xdo,640.*xdo.^.34,'-','color',rgb('indianred'),'linewidth',2);
% lstr=[lstr 'CNA; Erickson et al. (2004) [b=-0.5]'];
% 
% 
% %% for ENA
% 
% %Cramer 2012
% loglog(xdo,550.*xdo.^.6,'-','color',rgb('purple'),'linewidth',2);
% lstr=[lstr 'ENA; Cramer (2012) [b=-0.5]'];
% 
% %AB95
% loglog(xdo,680.*xdo.^.36,'-','color',rgb('forestgreen'),'linewidth',2);
% lstr=[lstr 'ENA; Atkinson and Boore (1995) [b=-0.5]'];
% 
% %AB14
% loglog(xdo,525.*xdo.^.45,'-','color',rgb('cyan'),'linewidth',2);
% lstr=[lstr 'ENA; Atkinson and Boore (2014) [b=-0.5]'];
% 
% %BS11
% loglog(xdo,410.*xdo.^.5,'-','color',rgb('magenta'),'linewidth',2);
% lstr=[lstr 'ENA; Boatwright and Seekins (2011) [b=-0.5]'];
% 
% %Atkinson 2004
% loglog(xdo,max(repmat(1000,size(xdo)),893.*xdo.^.32),'-','color',rgb('orange'),'linewidth',2);
% lstr=[lstr 'ENA; Atkinson (2004) [b=-0.5]'];
% 
% %Erickson 2004
% loglog(xdo,650.*xdo.^.36,'-','color',rgb('brown'),'linewidth',2);
% lstr=[lstr 'ENA; Erickson et al. (2004) [b=-0.5]'];
% 
% 
% plot(99,9999,'wo','markerfacecolor','w'); lstr=[lstr newline];
% 
% %% for Appalachian
% 
% %mine
% Q=460;
% eta=0.51;
% loglog(xdo,Q.*xdo.^eta,':','color',rgb('black'),'linewidth',2.5);
% lstr=[lstr 'Appalachian; This Study [b=-0.5]'];
% 
% %Shi 96
% loglog(xdo,573.5.*xdo.^.465,':','color',rgb('dodgerblue'),'linewidth',2);
% lstr=[lstr 'Appalachian; Shi et al. (1996) [b=n/a]'];
% 
% 
% 
% %% finish the figure
% 
% xlabel('\itf\rm (Hz)')
% ylabel(['\itQ(f)'])
% axis([.8 20 100 yma])
% % title('Region 1: Gulf Coast')
% % text(.67,5e4,tsr,'fontsize',20);
%  
% lh=legend(lstr,'location','northwest'); set(lh,'fontsize',13)
% ah=gca; set(ah,'fontsize',13)
% 
% removewhitespace(gca,1);
% jeffsavefig_r300(fh1,['Q_Comparison.png']);
disp('finished')
    