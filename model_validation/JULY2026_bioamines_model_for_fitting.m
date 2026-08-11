close all;
clear;
clc

%% Read in Data
samera_data0 = readtable('samera_data.csv','ReadVariableNames',true);
samera_data0.falciparum_inc_per_1000(isnan(samera_data0.falciparum_inc_per_1000))=0;

% Case Data to an array
samera_data = samera_data0.falciparum_inc_per_1000;

% 6 Possible Parameter Estimates
% dd,  fraction of cases severe, recovery time, immunity time, fraction susceptible, # initial cases 
%par0 = [6.5 .5 100 30 .5 100];
%ubs = [8 0.5 500 500 1 5000];
%lbs = [1 0   0   0  .75   1];


%par0 =[6.8   0.2   100     50     3  1000];
par0 = [7.4 0.1 304 308 1 1500];

ubs = [9  0.2   1000    500     1   1500];
lbs = [6   0.2   0      0     1   1500];

%
%par0 = [8 .5];

%ubs = [10 1];
%lbs = [1 0];

t0_vec = 4;
param_vec = zeros(length(t0_vec),length(par0));
J_vec = zeros(length(t0_vec),1);
F_vec = J_vec;

out_struct = cell(length(t0_vec),1);
mm_struct = cell(length(t0_vec),1);

%%
for i=1:length(t0_vec)

    tq = tic;
    t0 = t0_vec(i);

    disp(t0);


    %[P,J,R,F] = lsqnonlin(@(pars) model_fitter(pars,samera_data, t0), par0, lbs, ubs);
    
    problem = createOptimProblem('lsqnonlin','x0',par0, 'lb',lbs,'ub',ubs,'objective',@(pars) model_fitter(pars,samera_data, t0));
    ms = MultiStart('UseParallel',true);
    [P,J,F,out,manymin] = run(ms, problem, 5);

    %[P,J,~,F] = lsqnonlin(@(pars) model_fitter(pars, samera_data, t0), par0, lbs, ubs);


    param_vec(i,:) = P;
    J_vec(i,:) = J;
    F_vec(i,:) = F;

    out_struct{i} = out;
    mm_struct{i} = manymin;

    tq2 = toc(tq);
    disp(['Update: ' num2str(tq2) ' seconds required for t_0 = ' num2str(t0)])

end


function dy = model_fitter(pars, samera_data, t0)

    model_output = fitting_function(pars);
    
    dy = 0.1*model_output(1:length(t0:43))-samera_data(t0:43);


end





%%
figure(11)
clf;
plot(1:43,samera_data(1:43),'-ok','markersize',4,'markerfacecolor','k')

hold on;

for j=1:length(t0_vec)

    cases_model = fitting_function(param_vec(j,:));

    plot((t0_vec(j)):43, 0.1*cases_model(1:length(t0_vec(j):43)),'linewidth',2)

end

hold off;

plot_dates = samera_data0.month;

set(gca,'xtick',2:4:44,'xticklabel',datestr(plot_dates(2:4:44),'dd-mmm '),'fontsize',14)

ylabel([{'reported {\itP. falciparum} malaria incidence';'Semera'}],'HorizontalAlignment','center')
xlabel('month')
xlim([4 44])


% The Below Function contains the same contents as
% baseline_biogenic_amine_model_explorer.m and includes a calculation of
% monthly incidence.
function incidence_out = fitting_function(pars)

    
    p.DAYS = 10*365;     
    p.AGES = 50;       
    
    p.HumanPop = 1.8e6;
    
    %% Scenario Matrix
    
   
    % Vary: betaVH / betaHV / fecundity / bf / inf bf /biting rate
    % choice_mat = [0 0 0 0 0 0];    % Nothing affected by H/S
    choice_mat = [1 1 1 1 1 1];    % Everything affected by H/S
        
    nL = size(choice_mat,1);
    
    p.a = 10^-pars(1);
    
    p.prob_susc = pars(5);

    sev_vec = pars(2);

    p.recovery_time = pars(3);
    
    p.IntroCase = pars(6);

    p.immunity_time = pars(4);   
    
    n_param = length(sev_vec); 
    
    % Vary params. 1-yes; 0-no
    % bf_vary = 1;     
    % f_vary = 1;
    % betaVH_vary = 0;
    % betaHV_vary = 0;
    % bf_vary_inf = 0;

    
    %% Fecundity Matrices
    fecundity_sim_params = {'Control','1nM H','10 nM H', '0.15 mM S', '1.5 mM S', ...
        '0.15 mM S / 10 nM H', '1.5 mM S / 1 nM H'};
    

  
    frac_fecund_mat = [.505 .772 .806;
                         .515 .742 .757;
                         .415 .765 .761;
                         .555 .737 .811;
                         .405 .673 .758;
                         .425 .80 .79;
                         .495 .766 .835]; 
    
    no_eggs_mat = [40.76237624	39.74736842	47.25925926;
                     38.99029126	48.53932584	53.47169811;
                     39.01204819	42.32673267	46.23529412;
                     39.93693694	48.57142857	46.48051948;
                     39.59259259	44.52777778	44.12765957;
                     39.10588235	46.66346154	46.328125;
                     39.93939394	49.13265306	46.37878788];
    
    biting_mat = [1 0.8248 0.6239 0.1410 0.1871; % healthy week 1-5
                    1.1851 1.1532 0.7901 0.4641 0.1355]; % malaria week 1-5
    
   
    blood_feeding4 = [1/2*(196/251+264/345) 214/280 273/331 289/391 191/280 281/345 210/268];

    blood_feeding14 = [(1/2)*(171/221+136/219) 166/221 136/218 139/208 154/225 122/206 177/243];

    blood_feeding4inf = [.69 .57 .43 .62 .58 .63 .57];
    blood_feeding11inf = [.74 .77 .72 .53 .80 .55 .70];
    
    
    % Control Data on Weekly Feeding
    % Weeks 1-7
    p.wklyfeed = [1 0.890207715 0.606217617 0.597315436 0.539215686 0.530612245 0.714285714];
    
    
    %% Params
    
    % Human Values Estimated from Literature
    % Pre-patent (time from infection to detectable parasitemia: 10 days (5-10)
    % Incubation period of P. fal (infection to development of symptoms: 11
    % (6-14)
    % Non p.val incubation periods are longer (15-16 days)
    p.sigmaH = 1-exp(-1/13);       % 10 days human Incubation
    
    
    % Bretscher et al. Estimates Plasmodium duration with exp function:
    % Seasonally varying FOI: 43.9, 20, 18.9, 4.5, 14.9, 40.7 days
    % Duration of Infection: 201 days? 
    % Filipe et al. 180 days
       
    p.gammaH = 1-exp(-1/p.recovery_time);       % 200 days recovery
    
    p.gammaH2 = 1-exp(-1/p.recovery_time);      % Severe malaria Recovery
    
    % Duration of Temporary Immunity
    % 300 days
    p.psiH = (1-exp(-1/p.immunity_time));
    
    p.psiH = (1-exp(-1/p.immunity_time));
    
    
    % Percentage of severe malaria cases
    p.sevMal = 0.1;
    
    
    % Vector Values Estimated from Literature
    % Filipe et al. Use 0.25 (prob inco on bite)
    % 
    %betaHV0 = .25;      % transmission from Human to Vector
    %betaHV02 = betaHV0;
    %betaHV2 = p.betaHV0;
    
    betaHV = 0.5652;  % Control
    betaHV02 = 0.6458; % Malaria
    
    %betaHV02 = 0.3721;
    betaHV0 = 0.3721; % Healthy
    
    
    % Filipe et al. use 0.3-0.4, 0.03, and 0.015 for Symptomatic, Asymp,
    % Undetecta
    p.betaVH0 = .4;      % transmission from Vector to Human
    
    betaVH0 = p.betaVH0;
    betaVH02 = p.betaVH0;
    
    
    % Biting Rate 
    % Filipe Et al. .6 per day
    % 2 is from Okuneye paper
    p.c = 3;
    
    % Density Dependent Parameters
    % p.a = .2e-6;
    %p.a = 1e-7;
    %p.a = 1e-6;
    p.b = 3.4;
    
    % Vector Population Parameters estimated from literature
    p.muE = .99;                % 99% of eggs survive to the next day
    p.hatch = 1-exp(-1/4.5);      % 4.5 days to hatch
    p.nuJV = 1-exp(-1/12);      % 12 days emergence time
    p.sigmaV = 1-exp(-1/10);    % 10 day EIP
    
    p.muJV = exp(-1/14);        % Daily Survival Probability
    
    
    % Age-Dependent Survival
    p.muA = [1	0.9838489536	0.9257992031	0.9495706861	0.8914462427	0.9495706861	...
           0.8914462427	0.932898906	0.8552311148	0.8756304787	0.810703546	0.8303462267	...
           0.810703546	0.8303462267	0.7740994629	0.8206801169	0.6853083172	0.8097616291 ...
   	    0.6425933253	0.8097616291	0.6425933253	0.7811960479	0.616879113	0.6637882967	...
        0.5850950644	0.5946303558	0.5850950644	0.5946303558	0.5487472736	0.573068736	...
        0.4162667513	0.5338247303	0.3442508055	0.5338247303	0.3442508055	0.5130757385 ...
        0.3379432651	0.3591177702	0.3018966395	0.2974613636	0.3018966395	0.2974613636 ...
        0.2804840208	0.2974613636	0.1738698697	0.2788483287	0.1006131485	0.2788483287 ...
        0.1006131485	0.2366442573	0.07042920397	0.1014189674	0.07042920397	0.06761264495 ...
        0.07042920397	0.06761264495	0.03521460199	0];
    
    
        %% Multi-Run Output Matrices
    
        
        prevalence = zeros(nL, p.DAYS, n_param);
        incidence = zeros(nL, p.DAYS, n_param);
        total_cases = zeros(nL, p.DAYS, n_param);
    
        prevalence2 = prevalence;
        incidence2 = incidence;
        total_cases2 = total_cases;
    
        EGGSSV = prevalence;
        EGGS01 = prevalence;
        EGGS02 = prevalence;
    
        infectedv1 = prevalence;
        infectedv2 = prevalence;
    
        NV1 = EGGSSV;
        NV2 = NV1;
    
        VHR = zeros(nL, p.DAYS, n_param);
        larval_pop = zeros(nL, p.DAYS, n_param);
    
        feeding = zeros(nL, p.AGES, n_param);
    
        p.fecundity_storage = zeros(nL,p.AGES);
        p.eggs_storage = zeros(nL,p.AGES);
    
    
    for ii=1:n_param
        p.sevMal = sev_vec(ii);
    
        for jj=1:nL
        
        
            betaVH_vary = choice_mat(jj,1);
            betaHV_vary = choice_mat(jj,2);
            f_vary = choice_mat(jj,3);
            bf_vary = choice_mat(jj,4);
            bf_vary_inf = choice_mat(jj,5);
    
            biting_vary = choice_mat(jj,6);
            
            p.cc = zeros(p.AGES,1);
            p.cc2 = p.cc;
    
            if biting_vary ==1
                p.cc(1:7) = p.c*biting_mat(1,1);
                p.cc(8:14) = p.c*biting_mat(1,2);
                p.cc(15:21) = p.c*biting_mat(1,3);
                p.cc(22:28) = p.c*biting_mat(1,4);
                p.cc(29:p.AGES) = p.c*biting_mat(1,5);
    
                p.cc2(1:7) = p.c*biting_mat(2,1);
                p.cc2(8:14) = p.c*biting_mat(2,2);
                p.cc2(15:21) = p.c*biting_mat(2,3);
                p.cc2(22:28) = p.c*biting_mat(2,4);
                p.cc2(29:p.AGES) = p.c*biting_mat(2,5);
            else
                p.cc(1:7) = p.c*biting_mat(1,1);
                p.cc(8:14) = p.c*biting_mat(1,2);
                p.cc(15:21) = p.c*biting_mat(1,3);
                p.cc(22:28) = p.c*biting_mat(1,4);
                p.cc(29:p.AGES) = p.c*biting_mat(1,5);
    
                p.cc2=p.cc;
            end
    
    
    
            %% Do beta values vary? 
            if betaVH_vary==1
                p.betaVH = betaVH0;
                p.betaVH2 = betaVH02;
            else
                p.betaVH = p.betaVH0;
                p.betaVH2 = p.betaVH0;
            end
        
            if betaHV_vary==1
                p.betaHV = betaHV0;
                p.betaHV2 = betaHV02;
            else
                p.betaHV = betaHV0;
                p.betaHV2 = betaHV0;
            end
        
        
            %% Fecundity Values
            %f_choice = 1;
            f_choice = 7;
            f_choice2 = 6;
            p.frac_fecund = zeros(1,p.AGES);
            p.no_eggs = zeros(1,p.AGES);
        
            p.frac_fecund2 = p.frac_fecund;
            p.no_eggs2 = p.no_eggs;
            
    
            %Modified for Shorter GC
            p.frac_fecund(2:6) = frac_fecund_mat(f_choice,1);
            p.frac_fecund(7:11) = frac_fecund_mat(f_choice,2);
            p.frac_fecund(12:p.AGES) = frac_fecund_mat(f_choice,3);
            
            p.no_eggs(2:6) = no_eggs_mat(f_choice,1)/2;
            p.no_eggs(7:11) = no_eggs_mat(f_choice,2)/2;
            p.no_eggs(12:p.AGES) = no_eggs_mat(f_choice,3)/2;
        
            
            if(f_vary==1)
    
                p.frac_fecund2(2:6) = frac_fecund_mat(f_choice2,1);
                p.frac_fecund2(7:11) = frac_fecund_mat(f_choice2,2);
                p.frac_fecund2(12:p.AGES) = frac_fecund_mat(f_choice2,3);
                
                p.no_eggs2(2:6) = no_eggs_mat(f_choice2,1)/2;
                p.no_eggs2(7:11) = no_eggs_mat(f_choice2,2)/2;
                p.no_eggs2(12:p.AGES) = no_eggs_mat(f_choice2,3)/2;
    
    
    
            else
                p.frac_fecund2 = p.frac_fecund;
                p.no_eggs2 = p.no_eggs;
            end
        
          
        
            %% Does Blood Feeding change? 
            % Blood Meal Vector
            p.av = zeros(p.AGES,1);
            p.av2 = p.av;
        
            % Second Blood Meal at Day 4
            p.av(1) = 0;

            dist = [.25 .5 .25];

            dist = [0.333 0.333 0.333];
        
            bf_num = 7;
            bf_num2 = 6;
        
            % first bloodmeal
            p.av(2:4) = blood_feeding4(bf_num)*dist;
        
            p.av(7:9) = 0.5*(blood_feeding4(bf_num)+blood_feeding14(bf_num))*dist;
        
            % second bloodmeal
            p.av(12:14) = blood_feeding14(bf_num)*dist;
        
            % 3rd and 4th Bloodmeals
            p.av(17:19) = .778*p.wklyfeed(3)*dist;
    
            % 4th and beyond
            p.av(22:24) = .778*p.wklyfeed(4)*dist;
            p.av(27:29) = .778*p.wklyfeed(4)*dist;
            p.av(32:34) = .778*p.wklyfeed(4)*dist;
            p.av(37:39) = .778*p.wklyfeed(4)*dist;
            p.av(42:44) = .778*p.wklyfeed(4)*dist;
            p.av(47:49) = .778*p.wklyfeed(4)*dist;
    
        
          if bf_vary==1
        
            % first 3 bloodmeals
            p.av2(2:4) = blood_feeding4(bf_num2)*dist;
            p.av2(7:9) = 0.5*(blood_feeding4(bf_num2)+blood_feeding14(bf_num2))*dist;
            p.av2(12:14) = blood_feeding14(bf_num2)*dist;
        
            % 5-6 Bloodmeals
            p.av2(17:19) = .778*p.wklyfeed(3)*dist;
            p.av2(22:24) = .778*p.wklyfeed(4)*dist;
            p.av2(27:29) = .778*p.wklyfeed(4)*dist;
            p.av2(32:34) = .778*p.wklyfeed(4)*dist;
            p.av2(37:39) = .778*p.wklyfeed(4)*dist;
            p.av2(42:44) = .778*p.wklyfeed(4)*dist;
            p.av2(47:49) = .778*p.wklyfeed(4)*dist;
           
          else 
              p.av2 = p.av;
        
          end
        
           %% Does feeding change with infection? 
           bfi_num = 7;
           bfi_num2 = 6;
           
           p.avi = p.av;
           p.avi2 = p.av;
           
           p.avi(2:4) = blood_feeding4inf(bfi_num)*dist;
           p.avi(7:9) = (0.5)*(blood_feeding4inf(bfi_num)+blood_feeding11inf(bfi_num))*dist;
           p.avi(12:14) = blood_feeding11inf(bfi_num)*dist;
           p.avi(17:19) = .778*p.wklyfeed(3)*dist;
           p.avi(22:24) = .778*p.wklyfeed(4)*dist;
           p.avi(27:29) = .778*p.wklyfeed(4)*dist;
           
           
           if bf_vary_inf==1
               p.avi2(2:4) = blood_feeding4inf(bfi_num2)*dist;
               p.avi2(7:9) = (0.5)*(blood_feeding4inf(bfi_num2)+blood_feeding11inf(bfi_num2))*dist;
               p.avi2(12:14) = blood_feeding11inf(bfi_num2)*dist;
               
               p.avi2(17:19) = .778*p.wklyfeed(3)*dist;
               p.avi2(22:24) = .778*p.wklyfeed(4)*dist;
               p.avi2(27:29) = .778*p.wklyfeed(4)*dist;
           else
                p.avi2 = p.avi;
           end
           
            
        
            %% Storage Setup
            
            % Human Compartments
            H = zeros(5, p.DAYS);
            
           
            % Make Decisions about Vectors habits
            % Day 0, Day 4, Day 14 Are Feeding Days
            % Day 1, 5, 15 are Resting Days
            % Day 2-3, 5-
            
            % Juvenile Vectors
            JV = zeros(2, p.DAYS);
            
            % Suceptible Vectors
            SV = zeros(p.AGES, p.DAYS);
            SV2 = SV;
            
            % Exposed Vectors
            EV = zeros(p.AGES, p.DAYS);
            EV2 = EV;
        
            % Infectious Vectors
            IV = zeros(p.AGES, p.DAYS);
            IV2 = IV;
        
            % Eggs Vector
            EGGS = zeros(1,p.DAYS);
            
            % Initialize
            NH0 = p.HumanPop;
            H(:,1) = [p.prob_susc*NH0 0 0 0 (1-p.prob_susc)*NH0];
            
            SV(1,1) = 100;
            JV(:,1) = 100;
            
           
            % Various Variables
            NH = zeros(1,p.DAYS);
            NV = zeros(1,p.DAYS);
            
            INFV = zeros(1,p.DAYS);
            lambdaVH = zeros(1,p.DAYS);
            lambdaHV11 = zeros(1,p.DAYS);
            lambdaHV12 = zeros(1,p.DAYS);
    
            INFV2 = zeros(1,p.DAYS);
            lambdaVH2 = zeros(1,p.DAYS);
            lambdaHV21 = zeros(1,p.DAYS);
            lambdaHV22 = zeros(1,p.DAYS);
    
            %% Run Model for Burn-in Period
            p.DAYS0=p.DAYS;
            
            output0 = discrete_eqns(p, NH, NV, H, SV, SV2, EV, IV, EV2, IV2, INFV, INFV2, lambdaHV11, lambdaHV12, lambdaHV21, lambdaHV22, lambdaVH, lambdaVH2, EGGS, JV);
    
    
            JV(:,1) = output0.JV(:,end);
            SV(:,1) = output0.SV(:,end);    
            
            H(:,1) = [p.prob_susc*NH0 0 0 0 (1-p.prob_susc)*NH0];
            H(3,1) = p.IntroCase;
            H(1,1) = H(1,1)-H(3,1);
            
            output = discrete_eqns(p, NH, NV, H, SV, SV2, EV, IV, EV2, IV2, INFV, INFV2, lambdaHV11, lambdaHV12, lambdaHV21, lambdaHV22, lambdaVH, lambdaVH2, EGGS, JV);
    
            
            %            % 
            prevalence(jj,:,ii) = output.H(3,:);
            incidence(jj,:,ii) = output.new_cases;
            total_cases(jj,:,ii) = cumsum(output.new_cases);
            
            prevalence2(jj,:,ii) = output.H(4,:);
            incidence2(jj,:,ii) = output.new_cases2;
            total_cases2(jj,:,ii) = cumsum(output.new_cases2);
        
            % 
            VHR(jj,:,ii) = output.NV./output.NH;
            % 
            larval_pop(jj,:,ii) = output.JV(2,:);
            % 
            feeding(jj,:,ii) = p.av;
    
    
            % Infected Proportion of Vectors
            infectedv1(jj,:,ii) = sum(output.EV+output.IV,1)./output.NV;
            infectedv2(jj,:,ii) = sum(output.EV2+output.IV2,1)./output.NV;
    
            % Look at Eggs
            EGGSSV(jj,:,ii) = sum(p.frac_fecund(:).*p.no_eggs(:).*(output.SV));
            EGGS01(jj,:,ii) = sum(p.frac_fecund(:).*p.no_eggs(:).*(output.EV+output.IV));
            EGGS02(jj,:,ii) = sum(p.frac_fecund2(:).*p.no_eggs2(:).*(output.EV2+output.IV2));
    
        
           % Total Vectors
           NV1(jj,:,ii) = sum(output.SV)+sum(output.EV)+sum(output.IV);
           NV2(jj,:,ii) = sum(output.SV2)+sum(output.EV2)+sum(output.IV2);
    
            
    
        end
    
    
    
    
    
    end
    
    %% Calculate Total Cases per Month and fraction of severe case
    y_0 = [31 28 31 30 31 30 31 31 30 31 30 31];
    y_L = [31 29 31 30 31 30 31 31 30 31 30 31];
    %L = [0 0 0 1 0 0 0 1 2];  % 0 for no leap year; 1 for leap year; 2 for final year which has only 7 months
    L = [0 0 1 2];
    yr = repmat(2017:2025,12,1);
    
    
    monthly_incidence = zeros(44,1);
    
    total_incidence = incidence+incidence2;
    
    t0 =1;
    
    for i=1:44  % i = 1, 13, 25,... is the first month of each year
    
        % this gives us the month associated with the month in i
        k = rem(i-1,12)+1;
    
        year_date = yr(i);
    
        if(rem(year_date,4)==0)
            days_in_month = y_L(k);
        else
            days_in_month = y_0(k);
        end
    
        temp = total_incidence(1,(t0:(t0+days_in_month)-1));
    
        monthly_incidence(i) = sum(temp);
    
    
        %display([i k year_date days_in_month t0 t0+days_in_month-1])
        
        t0 = t0+days_in_month;
    
    
        
    
    end
    
    incidence_out = monthly_incidence./p.HumanPop*1000;

end



%%
function output = discrete_eqns(p, NH, NV, H, SV, SV2, EV, IV, EV2, IV2, INFV, INFV2, lambdaHV11, lambdaHV12, lambdaHV21, lambdaHV22, lambdaVH, lambdaVH2, EGGS, JV)

    new_cases = zeros(1,p.DAYS);
    new_cases2 = new_cases;
    new_cases(1) = 0;

    for t=1:(p.DAYS-1)
    
        NH(t) = sum(H(:,t));
        NV(t) = sum(SV(:,t)+EV(:,t)+IV(:,t) + SV2(:,t) + IV2(:,t) + EV2(:,t));
        
        % Malaria
        INFV(t) = sum(p.cc.*p.avi.*IV(:,t));  
        lambdaVH(t) = 1-exp(-INFV(t)/NH(t));

        % Severe Malaria
        INFV2(t) = sum(p.cc2.*p.avi2.*IV2(:,t));        
        lambdaVH2(t) = 1-exp(-INFV2(t)/NH(t));
    
        % Human Dynamics
        H(1,t+1) = H(1,t)*(1-p.betaVH*lambdaVH(t)-p.betaVH2*lambdaVH2(t)) + p.psiH*H(5,t) ;
        H(2,t+1) = H(2,t)*(1-p.sigmaH) + (p.betaVH*lambdaVH(t)+p.betaVH2*lambdaVH2(t))*H(1,t);

        H(3,t+1) = H(3,t)*(1-p.gammaH) + (1-p.sevMal)*p.sigmaH*H(2,t);
        H(4,t+1) = H(4,t)*(1-p.gammaH2) + (p.sevMal)*p.sigmaH*H(2,t);

        H(5,t+1) = H(5,t)*(1-p.psiH) + p.gammaH*H(3,t) + p.gammaH2*H(4,t);
        
        new_cases(t+1) = (1-p.sevMal)*p.sigmaH*H(2,t);
        new_cases2(t+1) = p.sevMal*p.sigmaH*H(2,t);
       
    
    for a=1:(p.AGES-1)
    
        %lambdaHV(a,t) = 1-exp(-p.betaHV*p.c*p.av(a)*H(3,t)/NH(t));
        %lambdaHV(t) = 1-exp(-p.betaHV*p.cc(a)*H(3,t)/NH(t));


        % Type 1 vector bite a type 1 host
        lambdaHV11(t) = 1-exp(-p.cc(a)*H(3,t)/NH(t));
        
        % Type 1 vector bite a type 2 bite a type 2 host 
        lambdaHV12(t) = 1-exp(-p.cc(a)*H(4,t)/NH(t));

        % Type 2 vector bite a type 1 host
        lambdaHV21(t) = 1-exp(-p.cc2(a)*H(3,t)/NH(t));
        
        % Type 2 vector bite a type 2 bite a type 2 host 
        lambdaHV22(t) = 1-exp(-p.cc2(a)*H(4,t)/NH(t));

        % lambdaHV2(t) = 1-exp(-p.betaHV2*p.cc2(a)*H(4,t)/NH(t));
   
        % Healthy Level Bioamines
        %SV(a+1, t+1) = p.muA(a)*SV(a,t) - p.muA(a)*p.av(a)*(lambdaHV(t)+lambdaHV2(t))*SV(a,t);

        SV(a+1, t+1) = p.muA(a)*SV(a,t) - p.muA(a)*p.av(a)*SV(a,t)*(p.betaHV*lambdaHV11(t) + lambdaHV12(t));
            
        EV(a+1, t+1) = p.muA(a)*p.av(a)*lambdaHV11(t)*p.betaHV*SV(a,t) - p.muA(a)*p.sigmaV*EV(a,t) + p.muA(a)*EV(a,t);
    
        IV(a+1, t+1) = p.muA(a)*p.sigmaV*EV(a,t) + p.muA(a)*IV(a,t);

        % Severe Malaria Associated BioAmines
        %SV2(a+1, t+1) =  p.muA(a)*SV2(a,t) - p.muA(a)*p.av2(a)*(lambdaHV(t)+lambdaHV2(t))*SV2(a,t);
        SV2(a+1, t+1) =  p.muA(a)*SV2(a,t) + p.muA(a)*p.av(a)*lambdaHV12(t)*(1-p.betaHV2)*SV(a,t) - ...
                p.muA(a)*p.av2(a)*p.betaHV2*(lambdaHV21(t)+lambdaHV22(t))*SV2(a,t);


        EV2(a+1, t+1) = p.muA(a)*p.av(a)*lambdaHV12(t)*(p.betaHV2)*SV(a,t) + ...
                p.muA(a)*p.av2(a)*p.betaHV2*(lambdaHV21(t)+lambdaHV22(t))*SV2(a,t) - p.muA(a)*p.sigmaV*EV2(a,t) + p.muA(a)*EV2(a,t);
    
        IV2(a+1, t+1) = p.muA(a)*p.sigmaV*EV2(a,t) + p.muA(a)*IV2(a,t);
    
        EGGS(t) = EGGS(t)+p.frac_fecund(a)*p.no_eggs(a)*(SV(a,t) + EV(a,t) + IV(a,t)) ...
            + p.frac_fecund2(a)*p.no_eggs2(a)*(SV2(a,t) + EV2(a,t) + IV2(a,t));
    
    end
    
    
        % Juvenile and Newly Emerged Vector Dynamics
        
        JV(1,t+1) = p.muE*JV(1,t)*(1-p.hatch) + EGGS(t);
        JV(2,t+1) = p.muE*p.hatch*JV(1,t) + p.muJV*exp(-((p.a*JV(2,t))^p.b))*JV(2,t)*(1 - p.nuJV);
        
        SV(1, t+1) = (1/2)*p.nuJV*p.muJV*exp(-((p.a*JV(2,t))^p.b))*JV(2,t);
    
    end


    output.NH = NH;
    output.NV = NV;
    output.H = H;
    output.SV = SV;
    output.EV = EV;
    output.IV = IV;
    output.SV2 = SV2;
    output.EV2 = EV2;
    output.IV2 = IV2;
    output.INFV = INFV;
    output.lambdaVH = lambdaVH;
    output.lambdaHV11 = lambdaHV11;
    output.lambdaHV12 = lambdaHV12;
    output.lambdaVH2 = lambdaVH2;
    output.lambdaHV21 = lambdaHV21;
    output.lambdaHV22 = lambdaHV22;
    output.EGGS = EGGS;
    output.JV = JV;
    output.new_cases = new_cases;
    output.new_cases2 = new_cases2;

end