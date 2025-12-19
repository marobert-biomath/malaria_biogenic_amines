% Age Structured Malaria Model to Investigate Influences of Biogenic Amines
% on Malaria Transmission
% Author: Michael A. Robert, Virginia Tech
% Last Modified: 02 December 2025
close all;
clear 
clc

% Set Basic Parameters 
p.DAYS = 5*365;    % Maximum timespan of model simulations
p.AGES = 50;       % Maximum Age of Mosquitoes

p.IntroCase = 1;        % How many initial malaria cases will be introduced? 
p.HumanPop = 100000;    % What is the total human population size? 


%% Define a filename to figure output

%fileName = '003_feeding_only';
fileName = 'null';

%% Scenario Matrix

% Vary: [betaVH betaHV fecundity bloodfeeding_tendency infected_bft biting_rate
% The choice_mat vector is set up to run two simulations, one with the
% options designated in the first row, and a second with options designated
% in the second row. For each, 0 means the effect of biogenic amines has been turned off and 1
% means the effect has been turned on. 
% All zeros indicate that there are no impacts of malaria-associated biogenic 
% amines on any traits
% All ones indicate that all traits are influenced by malaria-associated biogenic amines
% Column 1: Is betaHV, transmission from hosts to vectors
% Column 2: betaVH, transmission from vectors to hosts
% Column 3: fecundity: both # of offspring per female and the % of females
%           laying eggs
% Column 4: bloodfeeding_tendency: the tendency to bloodfeed on each day of
%           uninfected mosquitoes
% Column 5: infected_bft: the tendenc to bloodfeed on each day of infected
%           mosquitoes
% Column 6: biting_rate: the biting rate as approximated by cirdadian
%           activity
%
choice_mat = [0 0 0 0 0 0; % No effects of H/S
              1 1 1 1 1 1; ...   % Everything affected by H/S
                        ];
nL = size(choice_mat,1);

bf_scen = {'No Effects';'With Effects'};

sev_vec = 0.2:.2:1;
n_param = length(sev_vec);

% Vary params. 1-yes; 0-no
% bf_vary = 1;     
% f_vary = 1;
% betaVH_vary = 0;
% betaHV_vary = 0;
% bf_vary_inf = 0;

%scen_name = 'combo_bf_fecundity_a_1e-7_';
%scen_name = 'fecundity_a_1e-7';
%scen_name = 'beta0p4-ctrl_mal_bf_a_1e-7';
%scen_name = 'multi_type_run_a_1e-7';

%% Fecundity Matrices
fecundity_sim_params = {'Control','1nM H','10 nM H', '0.15 mM S', '1.5 mM S', ...
    '0.15 mM S / 10 nM H', '1.5 mM S / 1 nM H'};

% xCols = [0 0 0;
%                     1 0 0;
%                     0.5 0.023 0.01;
%                     0 1 1;
%                     0 0 1;
%                     0.8 0.4 1;
%                     0.4 0 0.5];

xCols = lines;

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
biting_mat_inf_mult = [0.6762    1.0256    0.4826    0.4826    0.4826;
                0.7657    1.0957    0.5567    0.5567    0.5567];

biting_mat_inf = biting_mat.*biting_mat_inf_mult;

%bf_scen = {'Control','1nM H','10 nM H', '0.15 mM S', '1.5 mM S', ...
 %   '0.15 mM S / 10 nM H (Malaria)', '1.5 mM S / 1 nM H (Healthy)'};
%bf_scen = {'Control','1nM H','10 nM H', '0.15 \muM 5-HT', '1.5 \muM 5-HT', ...
%    'Severe Malaria Treatment', 'Healthy Treatment'};

blood_feeding4 = [1/2*(196/251+264/345) 214/280 273/331 289/391 191/280 281/345 210/268];

MOE4 = 1.96*sqrt(blood_feeding4.*(1-blood_feeding4)./[1/2*(251+345) 280 331 391 280 345 268]);


blood_feeding14 = [(1/2)*(171/221+136/219) 166/221 136/218 139/208 154/225 122/206 177/243];

MOE14 = 1.96*sqrt(blood_feeding14.*(1-blood_feeding14)./[1/2*(221+219) 221 218 208 225 206 243]);


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
p.gammaH = 1-exp(-1/200);       % 200 days recovery

p.gammaH2 = 1-exp(-1/200);      % Severe malaria Recovery

% Duration of Temporary Immunity
% 300 days
p.psiH = (1-exp(-1/300));

% Percentage of severe malaria cases
p.sevMal = 0.1;


% Vector Values Estimated from Literature
% 
%
% P. yoelii Experiments
betaHV = 0.5652;    % Control
betaHV02 = 0.6458;  % MABAC
betaHV0 = 0.3721;   % HABAC

% P. falciparium Experiments
betaHVf = 0.46;     % Control
betaHV02f = 0.6863; % MABAC 
betaHV0f = 0.4286;  % HABAC

% Transmission from Vector to Human 
% This is currently set up to be for all scenarios. 
p.betaVH0 = .4;      
betaVH0 = p.betaVH0;
betaVH02 = p.betaVH0;

% Biting Rate 
% Baseline biting rate
p.c = 3;

% Density Dependent Parameters
p.a = 1e-6;
p.b = 3.4;

% Vector Population Parameters estimated from literature
p.muE = .99;                % 99% of eggs survive each day
p.hatch = 1-exp(-1/4.5);    % 4.5 days: average time to hatch
p.nuJV = 1-exp(-1/12);      % 12 days: average time to emergence from egg to adult
p.sigmaV = 1-exp(-1/10);    % 10 days: average EIP
p.muJV = exp(-1/14);        % 14 days: average density-independent juvenile mosquito lifespan

% Age-Dependent Survival
% Taken from the control group of Coles et al. 2023
p.muA = [1	0.9838489536	0.9257992031	0.9495706861	0.8914462427	0.9495706861	...
        0.8914462427	    0.932898906	    0.8552311148	0.8756304787	0.810703546	    0.8303462267	...
        0.810703546	        0.8303462267	0.7740994629	0.8206801169	0.6853083172	0.8097616291 ...
   	    0.6425933253	    0.8097616291	0.6425933253	0.7811960479	0.616879113	    0.6637882967	...
        0.5850950644	    0.5946303558	0.5850950644	0.5946303558	0.5487472736	0.573068736	...
        0.4162667513	    0.5338247303	0.3442508055	0.5338247303	0.3442508055	0.5130757385 ...
        0.3379432651	    0.3591177702	0.3018966395	0.2974613636	0.3018966395	0.2974613636 ...
        0.2804840208	    0.2974613636	0.1738698697	0.2788483287	0.1006131485	0.2788483287 ...
        0.1006131485	    0.2366442573	0.07042920397	0.1014189674	0.07042920397	0.06761264495 ...
        0.07042920397	    0.06761264495	0.03521460199	0];


%% Set up matrices to store output from multiple runs
% nL is the number of rows in choice_mat, or the number of different
% simulations to consider
% n_param is the number of different values of the parameter that will be
% varied

% Stable Age Distribution of the population 
% and adult population
%SAD_output = zeros(nL, p.AGES+2);
%SAD_adult_output = zeros(nL, p.AGES);

% Prevalence, incidence, and total cases of non-severe malaria in humans
prevalence = zeros(nL, p.DAYS, n_param);
incidence = prevalence;
total_cases = prevalence;

% Prevalence, incidence, and total cases of severe malaria in humans
prevalence2 = prevalence;
incidence2 = prevalence;
total_cases2 = prevalence;

% Egg counts
% V - total number of eggs produced by S vectors
% 01 - total number of eggs produced by type HABAC vectors
% 02 - total number of eggs produced by type MABAC vectors
%EGGSSV = prevalence;
%EGGS01 = prevalence;
%EGGS02 = prevalence;

%infectedv1 = prevalence;
%infectedv2 = prevalence;

%NV1 = EGGSSV;
%NV2 = NV1;

%VHR = prevalence;
%larval_pop = prevalence;

%feeding = zeros(nL, p.AGES, n_param);

%p.fecundity_storage = zeros(nL,p.AGES);
%p.eggs_storage = zeros(nL,p.AGES);


for ii=1:n_param
    p.sevMal = sev_vec(ii);

    for jj=1:nL
    
        % Select trait values based on which scenario you want to run. 
        betaVH_vary = choice_mat(jj,1);
        betaHV_vary = choice_mat(jj,2);
        f_vary = choice_mat(jj,3);
        bf_vary = choice_mat(jj,4);
        bf_vary_inf = choice_mat(jj,5);
        biting_vary = choice_mat(jj,6);
        
        p.cc = zeros(p.AGES,1);
        p.cc2 = p.cc;
        p.cci = p.cc;
        p.cci2 = p.cci;

        % Set biting rate depending on scenario 
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

            p.cci(1:7) = p.c*biting_mat_inf(1,1);
            p.cci(8:14) = p.c*biting_mat_inf(1,2);
            p.cci(15:21) = p.c*biting_mat_inf(1,3);
            p.cci(22:28) = p.c*biting_mat_inf(1,4);
            p.cci(29:p.AGES) = p.c*biting_mat_inf(1,5);
            
            p.cc2i(1:7) = p.c*biting_mat_inf(2,1);
            p.cc2i(8:14) = p.c*biting_mat_inf(2,2);
            p.cc2i(15:21) = p.c*biting_mat_inf(2,3);
            p.cc2i(22:28) = p.c*biting_mat_inf(2,4);
            p.cc2i(29:p.AGES) = p.c*biting_mat_inf(2,5);

        else

            p.cc(1:7) = p.c*biting_mat(1,1);
            p.cc(8:14) = p.c*biting_mat(1,2);
            p.cc(15:21) = p.c*biting_mat(1,3);
            p.cc(22:28) = p.c*biting_mat(1,4);
            p.cc(29:p.AGES) = p.c*biting_mat(1,5);

            p.cci(1:7) = p.c*biting_mat_inf(1,1);
            p.cci(8:14) = p.c*biting_mat_inf(1,2);
            p.cci(15:21) = p.c*biting_mat_inf(1,3);
            p.cci(22:28) = p.c*biting_mat_inf(1,4);
            p.cci(29:p.AGES) = p.c*biting_mat_inf(1,5);

            p.cc2=p.cc;
            p.cc2i=p.cci;

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
        
        % p.frac_fecund(2:16) = frac_fecund_mat(f_choice,1);
        % p.frac_fecund(17:31) = frac_fecund_mat(f_choice,2);
        % p.frac_fecund(32:p.AGES) = frac_fecund_mat(f_choice,3);
        % 
        % p.no_eggs(2:16) = no_eggs_mat(f_choice,1)/2;
        % p.no_eggs(17:31) = no_eggs_mat(f_choice,2)/2;
        % p.no_eggs(32:p.AGES) = no_eggs_mat(f_choice,3)/2;

        %Modified for Shorter GC
        p.frac_fecund(2:6) = frac_fecund_mat(f_choice,1);
        p.frac_fecund(7:11) = frac_fecund_mat(f_choice,2);
        p.frac_fecund(12:p.AGES) = frac_fecund_mat(f_choice,3);
        
        p.no_eggs(2:6) = no_eggs_mat(f_choice,1)/2;
        p.no_eggs(7:11) = no_eggs_mat(f_choice,2)/2;
        p.no_eggs(12:p.AGES) = no_eggs_mat(f_choice,3)/2;
    
        
        if(f_vary==1)
            % p.frac_fecund2(2:16) = frac_fecund_mat(f_choice2,1);
            % p.frac_fecund2(17:31) = frac_fecund_mat(f_choice2,2);
            % p.frac_fecund2(32:p.AGES) = frac_fecund_mat(f_choice2,3);
            % 
            % p.no_eggs2(2:16) = no_eggs_mat(f_choice2,1)/2;
            % p.no_eggs2(17:31) = no_eggs_mat(f_choice2,2)/2;
            % p.no_eggs2(32:p.AGES) = no_eggs_mat(f_choice2,3)/2;

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
        %dist = [.05 .25 .4 .25 .05];
        dist = [.25 .5 .25];
        %dist = [0 1 0];
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
        
        % Juvenile Vectors
        JVEC = zeros(1, p.DAYS);
        
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
        H(:,1) = [NH0 0 0 0 0];
        
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
        
        H(3,1) = p.IntroCase;
        H(1,1) = H(1,1)-H(3,1);
      
        
        display([jj ii])
        
        % Run Model 
       
        output = discrete_eqns(p, NH, NV, H, SV, SV2, EV, IV, EV2, IV2, INFV, INFV2, lambdaHV11, lambdaHV12, lambdaHV21, lambdaHV22, lambdaVH, lambdaVH2, EGGS, JV);


        % Store outputs for each simulation
        prevalence(jj,:,ii) = output.H(3,:);
        incidence(jj,:,ii) = output.new_cases;
        total_cases(jj,:,ii) = cumsum(output.new_cases);
        
        prevalence2(jj,:,ii) = output.H(4,:);
        incidence2(jj,:,ii) = output.new_cases2;
        total_cases2(jj,:,ii) = cumsum(output.new_cases2);
    
        %% Optional Output for looking at different measures
        %VHR(jj,:,ii) = output.NV./output.NH;
        % 
        %larval_pop(jj,:,ii) = output.JV(2,:);
        % 
        %feeding(jj,:,ii) = p.av;


        % Infected Proportion of Vectors
        %infectedv1(jj,:,ii) = sum(output.EV+output.IV,1)./output.NV;
        %infectedv2(jj,:,ii) = sum(output.EV2+output.IV2,1)./output.NV;

        % Look at Eggs
        %EGGSSV(jj,:,ii) = sum(p.frac_fecund(:).*p.no_eggs(:).*(output.SV));
        %EGGS01(jj,:,ii) = sum(p.frac_fecund(:).*p.no_eggs(:).*(output.EV+output.IV));
        %EGGS02(jj,:,ii) = sum(p.frac_fecund2(:).*p.no_eggs2(:).*(output.EV2+output.IV2));

    
       % Total Vectors
       %NV1(jj,:,ii) = sum(output.SV)+sum(output.EV)+sum(output.IV);
       %NV2(jj,:,ii) = sum(output.SV2)+sum(output.EV2)+sum(output.IV2);

        

    end
end



%% Calculate Peak and Timing of Peak of Incidence

tt = 0:(p.DAYS-1);

max_inc = zeros(size(incidence,1),size(incidence,3));
time_max = max_inc;
duration = time_max;

for j=1:size(incidence,3)

    for i=1:size(incidence,1)       
    
        tot_inc = incidence(i,:,j)+incidence2(i,:,j);

        max_inc(i,j) = max(tot_inc);
    
        if(max_inc(i,j)~=0)
            time_max(i,j) = tt(tot_inc==max_inc(i,j));
            tmin = min(tt(tot_inc>0));
            tmax = max(tt(tot_inc>0));
            duration(i,j) = tmax-tmin;
        else
            time_max(i,j) = 0;
            duration(i,j) = 0;
        end
     
    end

end

%% Save Output as a .mat file for figures 
save(fileName) 


%% Plot for Paper 
% Please note that for visibility, some plot bounds may need to be
% modified depending on which scenarios are being run. 
a = 150;    b = 3.6*365;

ToC = 4*365;

bc1 = [.3 .3 .3];
bc2 = [.4 .4 .9];

xLimz = [.2-.15 1.15];

f3 = figure(1);
clf;
set(f3, 'Color','White','Units','Inches','PaperUnits','Inches','Position', [1 4 14 8], ...
    'PaperSize',[14 8])

tcl = tiledlayout(2,3);

ax1 = nexttile([1 2]);
hold on;
box on;
grid on;

ax3 = nexttile;
hold on;
box on;

ax4 = nexttile;
hold on;
box on;

ax5 = nexttile;
hold on;
box on;

ax6 = nexttile;
hold on;
box on;

for jjj=1:length(sev_vec)

    % Plot Incidence
    h(jjj)=plot(ax1,0:(p.DAYS-1), squeeze(incidence(1,:,jjj)+incidence2(1,:,jjj)), 'LineWidth',3,'Color',xCols(jjj+2,:));
    plot(ax1,0:(p.DAYS-1), squeeze(incidence(2,:,jjj)+incidence2(2,:,jjj)), '--','LineWidth',3,'Color',xCols(jjj+2,:))



end

leghandls = h;

lz = legend(h,[repmat('\delta = ',length(sev_vec),1) num2str(sev_vec')],'Box','off');
set(lz,'position',[.15 .94 .3 .3],'FontSize',14,'orientation','horizontal')

ax1; hold off;
xlabel(ax1, 'day')
ylabel(ax1,{'daily incidence';'per 100,000 people'})
set(ax1,'fontsize',14)
xlim(ax1,[a b])
ylim(ax1,[0 600])

xx1 = get(ax1,'XLim');
yy1 = get(ax1,'YLim');

x1x = xx1(1)+.02*(xx1(2)-xx1(1));
y1y = yy1(1)+.9*(yy1(2)-yy1(1));

text(ax1,x1x,y1y,'(a)','FontSize',14)

%%%%%%% PANEL (B) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Calculate % Change in cases from HABAC to including MABAC effects
prctChangeC = (temp(2,:)-temp(1,:))./(temp(1,:))*100;
b1 = bar(ax3,sev_vec,prctChangeC,'facecolor','flat');

ylabel(ax3,{'% change in total cases within';'four years of introduction'})

xlabel(ax3,'fraction of cases of severe malaria (\delta)')

ylim(ax3,[-20 8])
xlim(ax3,xLimz)
b1(1).CData = bc1;

set(ax3,'fontsize',14)
ax3; hold off;


xx1 = get(ax3,'XLim');
yy1 = get(ax3,'YLim');x1x = xx1(1)+.05*(xx1(2)-xx1(1));
y1y = yy1(1)+.9*(yy1(2)-yy1(1));

text(ax3,x1x,y1y,'(b)','FontSize',14)


temp = squeeze(prevalence(:,end-1,:)+prevalence2(:,end-1,:));

%%%%%%% PANEL (C) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate % Change in size of outbreak peak from HABAC to including MABAC effects
prctChangeP = (max_inc(2,:)-max_inc(1,:))./(max_inc(1,:))*100;

b1 = bar(ax4,sev_vec,prctChangeP,'facecolor','flat');

ylabel(ax4,{'% change in size';'of outbreak peak'})
xlabel(ax4,'fraction of cases of severe malaria (\delta)')
ylim(ax4,[-20 8])
xlim(ax4,xLimz)
set(ax4,'fontsize',14)

%ylim([0.9*min(max_inc) 1.1*max(max_inc)])
b1(1).CData = bc1;
%b1(2).CData = bc2;

xx1 = get(ax4,'XLim');
yy1 = get(ax4,'YLim');

x1x = xx1(1)+.05*(xx1(2)-xx1(1));
y1y = yy1(1)+.9*(yy1(2)-yy1(1));

text(ax4,x1x,y1y,'(c)','FontSize',14)


%%%%%%% PANEL (D) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate % Change in timing of outbreak peak from HABAC to including MABAC effects
prctChangeT = (time_max(2,:)-time_max(1,:))./(time_max(1,:))*100;

b1 = bar(ax5,sev_vec,prctChangeT,'facecolor','flat');
%b1=bar(ax5,sev_vec, time_max,'facecolor','flat');
ylabel(ax5,{'% change in time';'of outbreak peak'})
xlabel(ax5,'fraction of cases of severe malaria (\delta)')
%set(gca,'fontsize',14,'xticklabel',bf_scen(scen_vec)','XTickLabelRotation',90)
% set(gca,'fontsize',14,'xticklabel',bf_scen,'XTickLabelRotation',90)
set(ax5,'fontsize',14)

%ylim([0.9*min(time_max) 1.1*max(time_max)])
ylim(ax5,[0 20])
xlim(ax5,xLimz)
b1(1).CData = bc1;
%b1(2).CData = bc2;

xx1 = get(ax5,'XLim');
yy1 = get(ax5,'YLim');

x1x = xx1(1)+.05*(xx1(2)-xx1(1));
y1y = yy1(1)+.9*(yy1(2)-yy1(1));

text(ax5,x1x,y1y,'(d)','FontSize',14)

temp = squeeze(total_cases(:,ToC,:)+total_cases2(:,ToC,:));





%%%%%%% PANEL (E) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate % Change in endemic prevalence from HABAC to including MABAC effects
prctChangeE = (temp(2,:)-temp(1,:))./(temp(1,:))*100;
b1 = bar(ax6,sev_vec,prctChangeE,'facecolor','flat');


ylabel(ax6,{'% change in';'endemic prevalence'})

xlabel(ax6,'fraction of cases of severe malaria (\delta)')
ylim(ax6,[-8 5])
xlim(ax6,xLimz)
b1(1).CData = bc1;
set(ax6,'fontsize',14)
ax6; hold off;

xx1 = get(ax6,'XLim');
yy1 = get(ax6,'YLim');

x1x = xx1(1)+.05*(xx1(2)-xx1(1));
y1y = yy1(1)+.9*(yy1(2)-yy1(1));

text(ax6,x1x,y1y,'(e)','FontSize',14)

ax3; hold off;
ax4; hold off;


saveas(gcf,['Modeling Paper/' fileName '.pdf'])



%% Plots to look more closely at dynamics and values

a = 100;        b = 3*365;

% Figure 2
f3 = figure(2);
clf;
set(f3, 'Color','White','Units','Inches','PaperUnits','Inches','Position', [1 4 14 8], ...
    'PaperSize',[14 8])
%
tcl = tiledlayout(2,4);
ax1 = nexttile;
hold on;
ax2 = nexttile;
hold on;
ax3 = nexttile;
hold on;
ax4 = nexttile;
axis off;

ax5 = nexttile;
ax6 = nexttile;
ax7 = nexttile;
ax8 = nexttile;

for jjj=1:length(sev_vec)


    % Plot Prevalence
    h(jjj) = plot(ax1,0:(p.DAYS-1), squeeze(prevalence(1,:,jjj)+prevalence2(1,:,jjj)),'linewidth',2,'Color',xCols(jjj+2,:));
    plot(ax1,0:(p.DAYS-1), squeeze(prevalence(2,:,jjj)+prevalence2(2,:,jjj)),'--','linewidth',2,'Color',xCols(jjj+2,:))

    % Plot Incidence
    plot(ax2,0:(p.DAYS-1), squeeze(incidence(1,:,jjj)+incidence2(1,:,jjj)), 'LineWidth',2,'Color',xCols(jjj+2,:))
    plot(ax2,0:(p.DAYS-1), squeeze(incidence(2,:,jjj)+incidence2(2,:,jjj)), '--','LineWidth',2,'Color',xCols(jjj+2,:))

    % Plot Total Cases
    plot(ax3,0:(p.DAYS-1), squeeze(total_cases(1,:,jjj)+total_cases2(1,:,jjj)), 'LineWidth',2,'Color',xCols(jjj+2,:))
    plot(ax3,0:(p.DAYS-1), squeeze(total_cases(2,:,jjj)+total_cases2(2,:,jjj)), '--','LineWidth',2,'Color',xCols(jjj+2,:))

  

end

xlabel(ax1, 'day')
ylabel(ax1,{'Prevalence';})
set(ax1,'fontsize',14)
xlim(ax1,[a b])
ylim(ax1,[0 6e4])


xlabel(ax2, 'day')
ylabel(ax2,{'Daily Incidence';})
set(ax2,'fontsize',14)
xlim(ax2,[a b])
ylim(ax2,[0 800])

xlabel(ax3, 'day')
ylabel(ax3,{'Cumulative Cases'})
set(ax3,'fontsize',14)
xlim(ax3,[a b])
ylim(ax3,[0 2e5])

%lz = legend(bf_scen','orientation','horizontal','Box','off');

leghandls = h;

lz = legend(h,[repmat('\delta = ',length(sev_vec),1) num2str(sev_vec')],'orientation','vertical','Box','off','NumColumns',2);
%lz.Layout.Tile = 'North';
set(lz,'position',[.7 .65 .3 .3],'FontSize',14)


%print(['presentation_' scen_name '_comparison_cases_relative.png'],'-dpng')

ax1; hold off;
ax2; hold off;
ax3; hold off;
ax4; hold off;


b1=bar(ax5,sev_vec, max_inc,'facecolor','flat');
ylabel(ax5,{'size of peak of outbreak'})
xlabel(ax5,'fraction of cases of severe malaria (\delta)')
%set(gca,'fontsize',14,'xticklabel',bf_scen(scen_vec)','XTickLabelRotation',90)
% set(gca,'fontsize',14,'xticklabel',bf_scen,'XTickLabelRotation',90)
ylim(ax5,[0 600])

%ylim([0.9*min(max_inc) 1.1*max(max_inc)])
b1(1).CData = [.9 .6 .6];
b1(2).CData = [.4 .4 .9];



b1=bar(ax6,sev_vec, time_max,'facecolor','flat');
ylabel(ax6,{'time of peak'})
xlabel(ax6,'fraction of cases of severe malaria (\delta)')
%set(gca,'fontsize',14,'xticklabel',bf_scen(scen_vec)','XTickLabelRotation',90)
% set(gca,'fontsize',14,'xticklabel',bf_scen,'XTickLabelRotation',90)

%ylim([0.9*min(time_max) 1.1*max(time_max)])
ylim(ax6,[0 900])
b1(1).CData = [.9 .6 .6];
b1(2).CData = [.4 .4 .9];

b1=bar(ax7,sev_vec, squeeze(total_cases(:,end,:)+total_cases2(:,end,:)),'facecolor','flat');
ylabel(ax7,{'total cases'})
%set(gca,'fontsize',14,'xticklabel',bf_scen(scen_vec)','XTickLabelRotation',90)
% set(gca,'fontsize',14,'xticklabel',bf_scen,'XTickLabelRotation',90)
xlabel(ax7,'fraction of cases of severe malaria (\delta)')
%ylim([0.9*min(total_cases(:,end)) 1.1*max(total_cases(:,end))])
ylim(ax7,[0 6e5])
b1(1).CData = [.9 .6 .6];
b1(2).CData = [.4 .4 .9];

b1=bar(ax8,sev_vec, squeeze(prevalence(:,end-1,:)+prevalence2(:,end-1,:)),'facecolor','flat');
ylabel(ax8,{'endemic prevalence'})
%set(gca,'fontsize',14,'xticklabel',bf_scen(scen_vec)','XTickLabelRotation',90)
% set(gca,'fontsize',14,'xticklabel',bf_scen,'XTickLabelRotation',90)
%ylim([0.9*min(prevalence(:,end-1)) 1.1*max(prevalence(:,end-1))])
xlabel(ax8,'fraction of cases of severe malaria (\delta)')
ylim(ax8,[0 3.5e4])
b1(1).CData = [.9 .6 .6];
b1(2).CData = [.4 .4 .9];


lz2=legend('without MABAC effects','with MABAC effects','box','off');
set(lz2,'Position',[.75 .55 .2 .1],'fontsize',14)

%% Model Equations 
function output = discrete_eqns(p, NH, NV, H, SV, SV2, EV, IV, EV2, IV2, INFV, INFV2, lambdaHV11, lambdaHV12, lambdaHV21, lambdaHV22, lambdaVH, lambdaVH2, EGGS, JV)

    new_cases = zeros(1,p.DAYS);
    new_cases2 = new_cases;
    new_cases(1) = 0;

    for t=1:(p.DAYS-1)
    
        NH(t) = sum(H(:,t));
        NV(t) = sum(SV(:,t)+EV(:,t)+IV(:,t) + SV2(:,t) + IV2(:,t) + EV2(:,t));
        
        % Malaria
        INFV(t) = sum(p.cci.*p.avi.*IV(:,t));  
        lambdaVH(t) = 1-exp(-INFV(t)/NH(t));

        % Severe Malaria
        INFV2(t) = sum(p.cc2i.*p.avi2.*IV2(:,t));        
        lambdaVH2(t) = 1-exp(-INFV2(t)/NH(t));
    
        % Human Dynamics
        H(1,t+1) = H(1,t)*(1-p.betaVH*lambdaVH(t)-p.betaVH2*lambdaVH2(t)) + p.psiH*H(5,t);
        H(2,t+1) = H(2,t)*(1-p.sigmaH) + (p.betaVH*lambdaVH(t)+p.betaVH2*lambdaVH2(t))*H(1,t);
        H(3,t+1) = H(3,t)*(1-p.gammaH) + (1-p.sevMal)*p.sigmaH*H(2,t);
        H(4,t+1) = H(4,t)*(1-p.gammaH2) + (p.sevMal)*p.sigmaH*H(2,t);
        H(5,t+1) = H(5,t)*(1-p.psiH) + p.gammaH*H(3,t) + p.gammaH2*H(4,t);
        
        % New Cases (non-severe)
        new_cases(t+1) = (1-p.sevMal)*p.sigmaH*H(2,t);
        % New Severe Cases 
        new_cases2(t+1) = p.sevMal*p.sigmaH*H(2,t);
       
    
    for a=1:(p.AGES-1)

        % Type 1 vector bite a type 1 host
        lambdaHV11(t) = 1-exp(-p.cc(a)*H(3,t)/NH(t));
        
        % Type 1 vector bite a type 2 host
        lambdaHV12(t) = 1-exp(-p.cc(a)*H(4,t)/NH(t));

        % Type 2 vector bite a type 1 host
        lambdaHV21(t) = 1-exp(-p.cc2(a)*H(3,t)/NH(t));
        
        % Type 2 vector bite a type 2 host 
        lambdaHV22(t) = 1-exp(-p.cc2(a)*H(4,t)/NH(t));

   
        % Healthy Level Biogenic amines (HABAC)
        SV(a+1, t+1) = p.muA(a)*SV(a,t) - p.muA(a)*p.av(a)*SV(a,t)*(p.betaHV*lambdaHV11(t) + lambdaHV12(t));
            
        EV(a+1, t+1) = p.muA(a)*p.av(a)*lambdaHV11(t)*p.betaHV*SV(a,t) - p.muA(a)*p.sigmaV*EV(a,t) + p.muA(a)*EV(a,t);
    
        IV(a+1, t+1) = p.muA(a)*p.sigmaV*EV(a,t) + p.muA(a)*IV(a,t);

        % Severe Malaria Associated BioAmines (MABAC)
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