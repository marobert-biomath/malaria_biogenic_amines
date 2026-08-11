%% This code is used to generate Figures S2-S3. 
sev_vec = 0.2:0.2:1; % Fraction of cases that are severe 
ToC = 1460; % Time at which "total cases" are counted following introduction


%% Sensitivity Data - Immunity
imm{1} = load('sens_immunity_50.mat','total_cases','total_cases2','prevalence','prevalence2','max_inc','time_max','incidence','incidence2');
imm{2} = load('sens_immunity_100.mat','total_cases','total_cases2','prevalence','prevalence2','max_inc','time_max','incidence','incidence2');
imm{3} = load('sens_immunity_150.mat','total_cases','total_cases2','prevalence','prevalence2','max_inc','time_max','incidence','incidence2');
imm{4} = load('sens_immunity_200.mat','total_cases','total_cases2','prevalence','prevalence2','max_inc','time_max','incidence','incidence2');
imm{5} = load('sens_immunity_250.mat','total_cases','total_cases2','prevalence','prevalence2','max_inc','time_max','incidence','incidence2');
imm{6} = load('sens_immunity_300.mat','total_cases','total_cases2','prevalence','prevalence2','max_inc','time_max','incidence','incidence2');
imm{7} = load('sens_immunity_350.mat','total_cases','total_cases2','prevalence','prevalence2','max_inc','time_max','incidence','incidence2');
imm{8} = load('sens_immunity_400.mat','total_cases','total_cases2','prevalence','prevalence2','max_inc','time_max','incidence','incidence2');
imm{9} = load('sens_immunity_450.mat','total_cases','total_cases2','prevalence','prevalence2','max_inc','time_max','incidence','incidence2');
imm{10} = load('sens_immunity_500.mat','total_cases','total_cases2','prevalence','prevalence2','max_inc','time_max','incidence','incidence2');


%% Sensitivity Data - Recovery 
rec{1} = load('sens_recovery_NS_50.mat','total_cases','total_cases2','prevalence','prevalence2','max_inc','time_max','incidence','incidence2');
rec{2} = load('sens_recovery_NS_100.mat','total_cases','total_cases2','prevalence','prevalence2','max_inc','time_max','incidence','incidence2');
rec{3} = load('sens_recovery_NS_150.mat','total_cases','total_cases2','prevalence','prevalence2','max_inc','time_max','incidence','incidence2');
rec{4} = load('sens_recovery_NS_200.mat','total_cases','total_cases2','prevalence','prevalence2','max_inc','time_max','incidence','incidence2');



%% Effects of duration of Immunity 
imm_vec = 50:50:500;

imm_prctChange_cases = zeros(length(imm_vec),length(sev_vec));
imm_prctChange_prev = zeros(length(imm_vec),length(sev_vec));
imm_prctChange_max = zeros(length(imm_vec),length(sev_vec));
imm_prctChange_tMax = zeros(length(imm_vec),length(sev_vec));

for j=1:length(imm_vec)
    
    temp = squeeze(imm{j}.total_cases(:,ToC,:)+imm{j}.total_cases2(:,ToC,:));
    imm_prctChange_cases(j,:) = (temp(2,:)-temp(1,:))./(temp(1,:))*100;

    temp = squeeze(imm{j}.prevalence(:,ToC,:)+imm{j}.prevalence2(:,ToC,:));
    imm_prctChange_prev(j,:) = (temp(2,:)-temp(1,:))./(temp(1,:))*100;

    imm_prctChange_max(j,:) = (imm{j}.max_inc(2,:)-imm{j}.max_inc(1,:))./(imm{j}.max_inc(1,:))*100;
    imm_prctChange_tMax(j,:) = (imm{j}.time_max(2,:)-imm{j}.time_max(1,:))./(imm{j}.time_max(1,:))*100;


end

%% Effects of duration of non-severe infectious period
rec_vec = 50:50:200;

rec_prctChange_cases = zeros(length(rec_vec),length(sev_vec));
rec_prctChange_prev = zeros(length(rec_vec),length(sev_vec));
rec_prctChange_max = zeros(length(rec_vec),length(sev_vec));
rec_prctChange_tMax = zeros(length(rec_vec),length(sev_vec));

for j=1:length(rec_vec)
    
    temp = squeeze(rec{j}.total_cases(:,ToC,:)+rec{j}.total_cases2(:,ToC,:));
    rec_prctChange_cases(j,:) = (temp(2,:)-temp(1,:))./(temp(1,:))*100;

    temp = squeeze(rec{j}.prevalence(:,ToC,:)+rec{j}.prevalence2(:,ToC,:));
    rec_prctChange_prev(j,:) = (temp(2,:)-temp(1,:))./(temp(1,:))*100;

    rec_prctChange_max(j,:) = (rec{j}.max_inc(2,:)-rec{j}.max_inc(1,:))./(rec{j}.max_inc(1,:))*100;
    rec_prctChange_tMax(j,:) = (rec{j}.time_max(2,:)-rec{j}.time_max(1,:))./(rec{j}.time_max(1,:))*100;


end


%%
f1 = figure(1);
set(f1,'Color','White')
subplot(2,2,1)
imagesc(sev_vec, imm_vec,imm_prctChange_cases)
set(gca,'YDir','normal')
xlabel('fraction of cases of severe malaria (\delta)')
ylabel('duration of immunity (days)')
colorbar
set(gca,'XTick',sev_vec,'XTickLabel',sev_vec)
title({'% change in total cases within';'four years of introduction'})
rectangle('Position',[0 275 1.1 50],'Linewidth',2,'EdgeColor','k')

subplot(2,2,2)
imagesc(sev_vec, imm_vec,imm_prctChange_prev)
set(gca,'YDir','normal')
xlabel('fraction of cases of severe malaria (\delta)')
ylabel('duration of immunity (days)')
colorbar
set(gca,'XTick',sev_vec,'XTickLabel',sev_vec)
title({'% change in endemic prevalence'})
rectangle('Position',[0 275 1.1 50],'Linewidth',2,'EdgeColor','k')

subplot(2,2,3)
imagesc(sev_vec, imm_vec,imm_prctChange_max)
set(gca,'YDir','normal')
xlabel('fraction of cases of severe malaria (\delta)')
ylabel('duration of immunity (days)')
colorbar
set(gca,'XTick',sev_vec,'XTickLabel',sev_vec)
title({'% change in size of outbreak peak'})
rectangle('Position',[0 275 1.1 50],'Linewidth',2,'EdgeColor','k')

subplot(2,2,4)
imagesc(sev_vec, imm_vec,imm_prctChange_tMax)
set(gca,'YDir','normal')
xlabel('fraction of cases of severe malaria (\delta)')
ylabel('duration of immunity (days)')
colorbar
set(gca,'XTick',sev_vec,'XTickLabel',sev_vec)
title({'% change in time of outbreak peak'})
rectangle('Position',[0 275 1.1 50],'Linewidth',2,'EdgeColor','k')


%%
rectangle_pos = [0 175 1.1 50];
f2 = figure(2);

set(f2,'Color','White')
subplot(2,2,1)
imagesc(sev_vec, rec_vec,rec_prctChange_cases)
set(gca,'YDir','normal')
xlabel('fraction of cases of severe malaria (\delta)')
ylabel({'duration of infectiousness';'non-severe malaria (days)'})
colorbar
set(gca,'XTick',sev_vec,'XTickLabel',sev_vec)
title({'% change in total cases within';'four years of introduction'})
rectangle('Position',rectangle_pos,'Linewidth',2,'EdgeColor','k')

subplot(2,2,2)
imagesc(sev_vec, rec_vec,rec_prctChange_prev)
set(gca,'YDir','normal')
xlabel('fraction of cases of severe malaria (\delta)')
ylabel({'duration of infectiousness';'non-severe malaria (days)'})
colorbar
set(gca,'XTick',sev_vec,'XTickLabel',sev_vec)
title({'% change in endemic prevalence'})
rectangle('Position',rectangle_pos,'Linewidth',2,'EdgeColor','k')

subplot(2,2,3)
imagesc(sev_vec, rec_vec,rec_prctChange_max)
set(gca,'YDir','normal')
xlabel('fraction of cases of severe malaria (\delta)')
ylabel({'duration of infectiousness';'non-severe malaria (days)'})
colorbar
set(gca,'XTick',sev_vec,'XTickLabel',sev_vec)
title({'% change in size of outbreak peak'})
rectangle('Position',rectangle_pos,'Linewidth',2,'EdgeColor','k')

subplot(2,2,4)
imagesc(sev_vec, rec_vec,rec_prctChange_tMax)
set(gca,'YDir','normal')
xlabel('fraction of cases of severe malaria (\delta)')
ylabel({'duration of infectiousness';'non-severe malaria (days)'})
colorbar
set(gca,'XTick',sev_vec,'XTickLabel',sev_vec)
title({'% change in time of outbreak peak'})
rectangle('Position',rectangle_pos,'Linewidth',2,'EdgeColor','k')

