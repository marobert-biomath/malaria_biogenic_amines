%% proporation calcs
stand_error = @(p, n) sqrt(p.*(1-p)./n);

%% Yoelii Sporozoites
% Control / Health / Malaria
ys_n_total = [46 43 48];
ys_n_pos = [26 16 31];
ys_frac = ys_n_pos./ys_n_total;
ys_se = stand_error(ys_frac, ys_n_total);

%% Falciparium Sporozoites
fs_n_total = [50 49 51];
fs_n_pos = [23 21 35];
fs_frac = fs_n_pos./fs_n_total;
fs_se = stand_error(fs_frac, fs_n_total);

%% Model output data
yoel = load('2025Aug01_yoelli_.mat');
falc = load('2025Aug01_falc_.mat');

%% Sporozoite data
spo_data = readtable('sporozoites.xlsx','readvariablenames',1);
spo_data2 = table2array(spo_data);

%% Figures
% Colors for yoelii and falciparum figures

cyoel = [.5 .5 .9];
cfalp = [.1 .9 .8];

sev_vec = yoel.sev_vec;
p.DAYS = yoel.p.DAYS;

% Figures Parameters 
small_fig_width = .18;
small_fig_height = .15;

small_fig_width_2 = .25;
small_fig_width3 = .22;

c1_x = .1;      r1_y = .8;
c2_x = .4;      r2_y = .45;
c3_x = .72;     r3_y = .15;

c2_x2 = .43;
c3_x2 = .76;


f1 = figure(1);
clf;
set(f1, 'Color','White','Units','Inches','PaperUnits','Inches','Position', [2 2 8 10], ...
    'PaperSize',[8 10])

%%%%%%% PANEL (d) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax1 = axes('Position',[c1_x r2_y small_fig_width3 small_fig_height]);

errorbar(ax1,.9:1:2.9, ys_frac, ys_se,'o','linewidth',1.5,'markerfacecolor',cyoel,'color',cyoel,'markersize',5)
hold on; 
errorbar(ax1,1.1:1:3.1, fs_frac, fs_se,'s','linewidth',1.5,'markerfacecolor',cfalp,'color',cfalp,'markersize',5)
hold off;
ylabel({'fraction positive';'for sporozoites'})
set(ax1,'xtick',1:3,'XtickLabel',{'Control', 'HABAC','MABAC'},...
    'fontsize',12,'XTickLabelRotation',30)
ylim([0.2 .8])
xlim([0.6 3.4])

xx1 = get(ax1,'XLim');
yy1 = get(ax1,'YLim');x1x = xx1(1)+.05*(xx1(2)-xx1(1));
y1y = yy1(1)+.9*(yy1(2)-yy1(1));

text(ax1,x1x,y1y,'(d)','FontSize',12)

%%%%%%% PANEL (e) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ax2 = axes('Position',[c2_x2 r2_y small_fig_width3 small_fig_height]);
bar(ax2,1, [ys_frac(3)./ys_frac(2)],'facecolor',cyoel)
hold on;
bar(ax2,3, [fs_frac(3)./fs_frac(2)],'facecolor',cfalp)
hold off;
set(ax2,'xtick',[1 3],'xticklabel',{'{\it P. yoelii}','{\it P. falciparum}'},'XTickLabelRotation',30)
ylabel({'increased transmission';'due to MABAC'})
set(ax2,'fontsize',12)

lz=legend('{\it P. yoelii}','{\it P. falciparum}','box','off','Orientation','horizontal');
set(lz,'fontsize',10, 'position',[.35 .65 .01 .01])

xx1 = get(ax2,'XLim');
yy1 = get(ax2,'YLim');x1x = xx1(1)+.05*(xx1(2)-xx1(1));
y1y = yy1(1)+.9*(yy1(2)-yy1(1));

text(ax2,x1x,y1y,'(e)','FontSize',12)

N = [50 49 51];

no_inf_c = sum(spo_data2(1:50,2)~=0);
no_inf_h = sum(spo_data2(1:50,3)~=0);
no_inf_m = sum(spo_data2(1:50,4)~=0);

inf_vec = [no_inf_c no_inf_h no_inf_m];
no_vec = N-inf_vec;

frac_vec =zeros(3,2);
frac_vec(:,1) = inf_vec'./N';
frac_vec(:,2) = no_vec'./N';

%%%%%%% PANEL (a) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ax2 = axes('Position',[c1_x r1_y small_fig_width small_fig_height]);

b1 = bar(ax2,[.7 1.7 2.7]+.15,inf_vec,'barwidth',.3,'FaceColor','flat');
hold on;
b2 = bar(ax2,[1 2 3]+.15,no_vec,'barwidth',.3,'facecolor','flat');
hold off;

b2.FaceColor = [.6 .6 .6];
b1.FaceColor = cfalp;

set(ax2,'xtick',1:3,'xticklabel',...
    {'Control', 'HABAC','MABAC'},...
    'XTickLabelRotation',30,'Fontsize',12)
ylim([0 50])
ylabel('number of mosquitoes')
xlim([0.6 3.4])


bars0 = zeros(2,2,3,3);
bars0(1,2,1,3) = 1;
bars0(1,2,2,3) = 2;
plot_stats(ax2,bars0)

xx1 = get(ax2,'XLim');
yy1 = get(ax2,'YLim');x1x = xx1(1)+.05*(xx1(2)-xx1(1));
y1y = yy1(1)+.9*(yy1(2)-yy1(1));

text(ax2,x1x,y1y,'(a)','FontSize',12)


%%%%%%% PANEL (b) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ax3 = axes('Position',[c2_x r1_y small_fig_width small_fig_height]);

b3=bar(ax3,frac_vec.*100,'stacked','Horizontal','off','BarWidth',.5,'FaceColor','flat');
set(ax3,'xticklabel',...
    {'Control', 'HABAC','MABAC'},...
    'fontsize',12)
ylabel('% of mosquitoes')
lz = legend('infected','uninfected');
set(lz,'Orientation','horizontal','box','off','position',[.35 .98 .01 .01],'fontsize',10)
xlim([0 3.4])

b3(1).CData = cfalp;
b3(2).CData = [.6 .6 .6];


c1 = [.8 .1 .1]; % red
c2 = [0 .1 1];
cols = [0 0 0; c2; c1];

xx1 = get(ax3,'XLim');
yy1 = get(ax3,'YLim');x1x = xx1(1)+.05*(xx1(2)-xx1(1));
y1y = yy1(1)+.9*(yy1(2)-yy1(1));

text(ax3,x1x,y1y,'(b)','FontSize',12)

%%%%%%% PANEL (c) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ax4a = axes('Position',[c3_x r1_y small_fig_width_2 .5*small_fig_height]);

pos = 3;
vs1=Violin({spo_data2(spo_data2(:,2)>0,2)},1,...
    'HalfViolin','full',...% left, full
    'QuartileStyle','none',... % boxplot, none
    'DataStyle', 'scatter',... % scatter, none
    'ShowNotches', false,...
    'ShowMean', false,...
    'ShowMedian', true,...
    'ViolinColor', {cols(1,:)},...
    'MarkerSize', 10,...
    'Orientation', 'vertical');

vs2=Violin({spo_data2(spo_data2(:,3)>0,3)},2,...
    'HalfViolin','full',...% left, full
    'QuartileStyle','none',... % boxplot, none
    'DataStyle', 'scatter',... % scatter, none
    'ShowNotches', false,...
    'ShowMean', false,...
    'ShowMedian', true,...
    'ViolinColor', {cols(2,:)},...
    'MarkerSize', 10,...
    'Orientation', 'vertical');

vs3=Violin({spo_data2(spo_data2(:,4)>0,4)},3,...
    'HalfViolin','full',...% left, full
    'QuartileStyle','none',... % boxplot, none
    'DataStyle', 'scatter',... % scatter, none
    'ShowNotches', false,...
    'ShowMean', false,...
    'ShowMedian', true,...
    'ViolinColor', {cols(3,:)},...
    'MarkerSize', 10,...
    'Orientation', 'vertical');
ylim([0 40])

set(ax4a,'xtick',1:3,'XtickLabel',{'Control', 'HABAC','MABAC'},...
    'fontsize',12,'XTickLabelRotation',30)

xlim([0.6 3.4])
yL = ylabel({'sporozoite count'});
set(yL,'Position',[.1 40])

ax4b = axes('Position',[c3_x r1_y+.6*(small_fig_height) small_fig_width_2 .5*small_fig_height]);

vs1=Violin({spo_data2(spo_data2(:,2)>0,2)},1,...
    'HalfViolin','full',...% left, full
    'QuartileStyle','none',... % boxplot, none
    'DataStyle', 'scatter',... % scatter, none
    'ShowNotches', false,...
    'ShowMean', false,...
    'ShowMedian', true,...
    'ViolinColor', {cols(1,:)},...
    'MarkerSize', 10,...
    'Orientation', 'vertical');

vs2=Violin({spo_data2(spo_data2(:,3)>0,3)},2,...
    'HalfViolin','full',...% left, full
    'QuartileStyle','none',... % boxplot, none
    'DataStyle', 'scatter',... % scatter, none
    'ShowNotches', false,...
    'ShowMean', false,...
    'ShowMedian', true,...
    'ViolinColor', {cols(2,:)},...
    'MarkerSize', 10,...
    'Orientation', 'vertical');

vs3=Violin({spo_data2(spo_data2(:,4)>0,4)},3,...
    'HalfViolin','full',...% left, full
    'QuartileStyle','none',... % boxplot, none
    'DataStyle', 'scatter',... % scatter, none
    'ShowNotches', false,...
    'ShowMean', false,...
    'ShowMedian', true,...
    'ViolinColor', {cols(3,:)},...
    'MarkerSize', 10,...
    'Orientation', 'vertical');
ylim([40 250])
xlim([0.6 3.4])
set(ax4b,'Xtick','','xcolor',[1 1 1])

xx1 = get(ax4b,'XLim');
yy1 = get(ax4b,'YLim');x1x = xx1(1)+.05*(xx1(2)-xx1(1));
y1y = yy1(1)+.9*(yy1(2)-yy1(1));

text(ax4b,x1x,y1y,'(c)','FontSize',12)

a = 150;    b = 3.6*365;

ToC = 4*365;

bc1 = [.5 .5 .9];

bc2 = [.1 .9 .8];


xLimz = [.2-.15 1.15];


%%%%%%% PANEL (f) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yoel.temp = squeeze(yoel.total_cases(:,ToC,:)+yoel.total_cases2(:,ToC,:));
falc.temp = squeeze(falc.total_cases(:,ToC,:)+falc.total_cases2(:,ToC,:));

ax5 = axes('Position',[c3_x2 r2_y small_fig_width3 small_fig_height]);


prctChangeC_y = (yoel.temp(2,:)-yoel.temp(1,:))./(yoel.temp(1,:))*100;
prctChangeC_f = (falc.temp(2,:)-falc.temp(1,:))./(falc.temp(1,:))*100;

data =[prctChangeC_y; prctChangeC_f];

b1 = bar(ax5,yoel.sev_vec,data,'facecolor','flat');


ylabel({'% change in total cases within';'four years of introduction'})
xlabel({'fraction of cases of';'severe malaria (\delta)'})

ylim([0 90])
xlim(xLimz)
b1(1).CData = bc1;
b1(2).CData = bc2;

set(ax5,'fontsize',12,'XTickLabelRotation',90)
hold off;


xx1 = get(ax5,'XLim');
yy1 = get(ax5,'YLim');x1x = xx1(1)+.05*(xx1(2)-xx1(1));
y1y = yy1(1)+.9*(yy1(2)-yy1(1));

text(gca,x1x,y1y,'(f)','FontSize',12)


%%%%%%% PANEL (g) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yoel.temp = squeeze(yoel.prevalence(:,end-1,:)+yoel.prevalence2(:,end-1,:));
falc.temp = squeeze(falc.prevalence(:,end-1,:)+falc.prevalence2(:,end-1,:));

prctChangeE_y = (yoel.temp(2,:)-yoel.temp(1,:))./(yoel.temp(1,:))*100;
prctChangeE_f = (falc.temp(2,:)-falc.temp(1,:))./(falc.temp(1,:))*100;

data = [prctChangeE_y' prctChangeE_f'];

ax5 = axes('Position',[c1_x r3_y small_fig_width3 small_fig_height]);

b1 = bar(ax5,yoel.sev_vec,data,'facecolor','flat');

ylabel({'% change in';'endemic prevalence'})

xlabel({'fraction of cases';'of severe malaria (\delta)'})
ylim(ax5,[0 30])
xlim(ax5,xLimz)
b1(1).CData = bc1;
b1(2).CData = bc2;
set(ax5,'fontsize',12,'XTickLabelRotation',90)
hold off;

xx1 = get(gca,'XLim');
yy1 = get(gca,'YLim');

x1x = xx1(1)+.05*(xx1(2)-xx1(1));
y1y = yy1(1)+.9*(yy1(2)-yy1(1));

text(gca,x1x,y1y,'(g)','FontSize',12)


%%%%%%% PANEL (h) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prctChangeP_y = (yoel.max_inc(2,:)-yoel.max_inc(1,:))./(yoel.max_inc(1,:))*100;
prctChangeP_f = (falc.max_inc(2,:)-falc.max_inc(1,:))./(falc.max_inc(1,:))*100;

data = [prctChangeP_y' prctChangeP_f'];


ax5 = axes('Position',[c2_x2 r3_y small_fig_width3 small_fig_height]);

b1 = bar(ax5,yoel.sev_vec,data,'facecolor','flat');
ylabel({'% change in size';'of outbreak peak'})
xlabel({'fraction of cases';'of severe malaria (\delta)'})

ylim([0 120])
xlim(xLimz)
set(ax5,'fontsize',12,'XTick',0.2:0.2:1,'XTickLabel',0.2:0.2:1,'xticklabelrotation',90)

b1(1).CData = bc1;
b1(2).CData = bc2;

xx1 = get(ax5,'XLim');
yy1 = get(ax5,'YLim');

x1x = xx1(1)+.05*(xx1(2)-xx1(1));
y1y = yy1(1)+.9*(yy1(2)-yy1(1));

text(x1x,y1y,'(h)','FontSize',12)


%%%%%%% PANEL (i) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prctChangeT_y = (yoel.time_max(2,:)-yoel.time_max(1,:))./(yoel.time_max(1,:))*100;
prctChangeT_f = (falc.time_max(2,:)-falc.time_max(1,:))./(falc.time_max(1,:))*100;

data = [prctChangeT_y' prctChangeT_f'];

ax5 = axes('Position',[c3_x2 r3_y small_fig_width3 small_fig_height]);

b1 = bar(ax5,yoel.sev_vec,data,'facecolor','flat');
%b1=bar(ax5,sev_vec, time_max,'facecolor','flat');

ylabel({'% change in time';'of outbreak peak'})
xlabel({'fraction of cases';'of severe malaria (\delta)'})
%set(gca,'fontsize',14,'xticklabel',bf_scen(scen_vec)','XTickLabelRotation',90)
% set(gca,'fontsize',14,'xticklabel',bf_scen,'XTickLabelRotation',90)
set(ax5,'fontsize',12,'XTickLabelRotation',90)

%ylim([0.9*min(time_max) 1.1*max(time_max)])
ylim(gca,[-50 20])
xlim(gca,xLimz)
b1(1).CData = bc1;
b1(2).CData = bc2;

xx1 = get(ax5,'XLim');
yy1 = get(ax5,'YLim');

x1x = xx1(1)+.05*(xx1(2)-xx1(1));
y1y = yy1(1)+.9*(yy1(2)-yy1(1));

text(ax5,x1x,y1y,'(i)','FontSize',12)

%saveas(gcf,['Modeling Paper/' '003_yoelli_falcip' '.pdf'])
saveas(gcf,['003_combined_figure_yoel_falcip' '.pdf'])


