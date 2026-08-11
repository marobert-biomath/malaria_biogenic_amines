% Age Structured Malaria Model to Investigate Influences of Biogenic Amines
% on Malaria Transmission
% Author: Michael A. Robert, Virginia Tech
% Last Modified: 02 December 2025
close all;
clear 
clc

% Set Basic Parameters that will not be modified for this study 
p.DAYS = 5*365;    % Maximum timespan of model simulations
p.AGES = 50;       % Maximum Age of Mosquitoes

p.IntroCase = 1;        % How many initial malaria cases will be introduced? 
p.HumanPop = 100000;    % What is the total human population size? 


%% Define a filename to figure output
% This 'fileName' shows up in the title for the model output and should be
% descriptive of the model simulation performed. Use filename 'null' if not
% interested in saving the results.  This will save a file that will be
% overwritten every time this filename is used. 

%fileName = 'v003_feeding_only_results';
fileName = 'null';

%% Scenario Matrix

% The vector choice_mat determines which traits will be varied to include biogenic amine impacts.
% choice_mat = [betaVH betaHV fecundity bloodfeeding_tendency infected_bft biting_rate]
% The choice_mat vector is set up to run two simulations, one with the
% options designated in the first row, and a second with options designated
% in the second row. For each, 0 means the effect of biogenic amines has been turned off and 1
% means the effect has been turned on. 
% All zeros indicate that there are no impacts of malaria-associated biogenic 
% amines on any traits
% All ones indicate that all traits are influenced by malaria-associated biogenic amines
% Column 1: Is betaHV, transmission from hosts to vectors
                % 0 P. yoelii values; no difference between HABAC and MABAC 
                % 1 P. yoelii values; difference between HABAC and MABAC
                % 2 P. falciparum values; no difference between HABAC and MABAC
                % 3 P. falciparum values; difference between HABAC and MABAC
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

% This vector determines the fraction of malaria cases that are severe.
% The model will run one simulation with each entry in sev_vec.
sev_vec = 0.2:.2:1;
n_param = length(sev_vec);

%% Data taken from Coles et al. 2023 and Ochwedo et al. 2025

% These are the treatment names associated with data presented in Coles et al. 2023 
% Coles, Taylor A., et al. "Ingested histamine and serotonin interact to alter Anopheles stephensi 
% feeding and flight behavior and infection with Plasmodium parasites." 
% Frontiers in physiology 14 (2023): 1247316.
treatment_names = {'Control', ...
                    '1nM H', ...
                    '10 nM H', ...
                    '0.15 mM S', ...
                    '1.5 mM S', ...
                    '0.15 mM S / 10 nM H', ...
                    '1.5 mM S / 1 nM H'};


% This matrix is the data from Coles et al. 2023 for the fraction of female
% mosquitoes that lay eggs.  Rows indicate experimental treatments, as outlined in 'treatment_names' above.  
% For this study, only rows 6 (MABAC) and 7 (HABAC) are used. 
% Columns indicate weeks of adult female mosquito life (week 1, week 2, week 3). 
frac_fecund_mat = [.505 .772 .806;
                     .515 .742 .757;
                     .415 .765 .761;
                     .555 .737 .811;
                     .405 .673 .758;
                     .425 .80 .79;
                     .495 .766 .835]; 

% This matrix is the data from Coles et al. 2023 for the average number of eggs per
% mosquitoes over two days.  Rows indicate experimental treatments, as outlined in 'treatment_names' above.  
% For this study, only rows 6 (MABAC) and 7 (HABAC) are used. 
% Columns indicate weeks of adult female mosquito life (week 1, week 2, week 3). 
no_eggs_mat = [40.76	39.75	47.26;
                 38.99	48.54	53.47;
                 39.01	42.33	46.24;
                 39.94	48.57	46.48;
                 39.59	44.53	44.13;
                 39.11	46.66	46.33;
                 39.94	49.13	46.38];

% This matrix calculates the fraction of uninfected adult female mosquitoes that
% bloodfeed at 4 days.  Each entry represents a treatment presented in
% Coles et al. 2023.
blood_feeding4 = [1/2*(196/251+264/345) 214/280 273/331 289/391 191/280 281/345 210/268];

% This matrix calculates the fraction of uninfected adult female mosquitoes that
% bloodfeed at 14 days.  Each entry represents a treatment presented in
% Coles et al. 2023.
blood_feeding14 = [(1/2)*(171/221+136/219) 166/221 136/218 139/208 154/225 122/206 177/243];


% This matrix calculates the fraction of infected adult female mosquitoes that
% bloodfeed at 4 days.  Each entry represents a treatment presented in
% Coles et al. 2023.
blood_feeding4inf = [.69 .57 .43 .62 .58 .63 .57];

% This matrix calculates the fraction of infected adult female mosquitoes that
% bloodfeed at 11 days.  Each entry represents a treatment presented in
% Coles et al. 2023.
blood_feeding11inf = [.74 .77 .72 .53 .80 .55 .70];


% Fraction of uninfected adult female mosquitoes who take a blood meal each
% week for 7 weeks. 
p.wklyfeed = [1 0.89 0.61 0.60 0.54 0.53 0.71];

% Fraction of mosquitoes who become infected upon biting an infectious host
% P. yoelii Experiments (Coles et al. 2023)
betaHV = 0.5652;    % Control
betaHV02 = 0.6458;  % MABAC
betaHV0 = 0.3721;   % HABAC

% Fraction of mosquitoes who become infected upon biting an infectious host
% P. falciparum Experiments (Figure 8 of the main text)
betaHVf = 0.46;     % Control
betaHV02f = 0.6863; % MABAC 
betaHV0f = 0.4286;  % HABAC


% Age-Dependent Survival
% Estimated from the control group shown in Coles et al. 2023
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


% The following two matrices are relative biting rates adapated from Ochwedo et al. 2025. 
% Relative biting rates are calculated for uninfected (biting_mat) based on  scaling 
% the weekly activity rates to that of the uninfected mosquitoes in week 1.
% Relative biting rates for infected (biting_mat_inf) mosquitoes are taken
% from calculating relative activity levels of infected mosquitoes compared
% to uninfected mosquitos (e.g. week 1 activity of HABAC infected mosquitoes
% compared to week 1 activity of HABAC uninfected mosquitoes).  This data
% is described in more detail in the text. 

% Ochwedo, Kevin O., et al. "Regulation of diel locomotor activity and retinal responses of 
% Anopheles stephensi by ingested histamine and serotonin is temperature-and infection-dependent."
% PLoS Pathogens 21.4 (2025): e1013139.

% uninfected biting rates
% columns are weeks 1-5; rows are HABAC (row 1) and MABAC (row 2)
biting_mat = [1 0.8248 0.6239 0.1410 0.1871; % healthy week 1-5
                1.1851 1.1532 0.7901 0.4641 0.1355]; % malaria week 1-5


% Infected biting rates
% columns are weeks 1-5; rows are HABAC (row 1) and MABAC (row 2)
biting_mat_inf_mult = [0.6762    1.0256    0.4826    0.4826    0.4826;
                0.7657    1.0957    0.5567    0.5567    0.5567];
biting_mat_inf = biting_mat.*biting_mat_inf_mult;


%% Model Parameters Estimated from literature
% Citations for all parameters are listed in Table 1 of the main text. 

% Human Values Estimated from Literature
% Pre-patent (time from infection to detectable parasitemia: 10 days (5-10)
% Incubation period of P. fal (infection to development of symptoms: 11
% (6-14)
% Non p.val incubation periods are longer (15-16 days)
p.sigmaH = 1-exp(-1/13);       % 13 days human Incubation


% Bretscher et al. Estimates Plasmodium duration with exp function:
% Seasonally varying FOI: 43.9, 20, 18.9, 4.5, 14.9, 40.7 days
% Duration of Infection: 201 days? 
% Filipe et al. 180 days
p.gammaH = 1-exp(-1/200);       % 200 days to recovery for non-severe malaria

p.gammaH2 = 1-exp(-1/200);      % 200 days to recovery for severe malaria 

% Duration of Temporary Immunity
p.psiH = (1-exp(-1/300));       % 300 days


% Transmission from Vector to Human 
% This is currently set up to be the same for all scenarios. 
p.betaVH0 = .4;      
betaVH0 = p.betaVH0;
betaVH02 = p.betaVH0;

% Biting Rate 
% Baseline biting rate - this is multiplied by the HABAC/MABAC biting rates
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



%% Set up matrices to store output from multiple runs
% nL is the number of rows in choice_mat, or the number of different
% simulations to consider
% n_param is the number of different values of the parameter that will be
% varied

% Prevalence, incidence, and total cases of non-severe malaria in humans
prevalence = zeros(nL, p.DAYS, n_param);
incidence = prevalence;
total_cases = prevalence;

% Prevalence, incidence, and total cases of severe malaria in humans
prevalence2 = prevalence;
incidence2 = prevalence;
total_cases2 = prevalence;


%% For loops to run through the scenarios with/without biogenic amine ffects (nL) and different
% fractions of severe malaria (n_param)

for ii=1:n_param
    p.sevMal = sev_vec(ii);

    for jj=1:nL
    
        % Select trait values jj based on which scenario you want to run. 
        betaVH_vary = choice_mat(jj,1);
        betaHV_vary = choice_mat(jj,2);
        f_vary = choice_mat(jj,3);
        bf_vary = choice_mat(jj,4);
        bf_vary_inf = choice_mat(jj,5);
        biting_vary = choice_mat(jj,6);
        
        % Set up biting rates depending on scenarios
        p.cc = zeros(p.AGES,1);
        p.cc2 = p.cc;
        p.cci = p.cc;
        p.cci2 = p.cci;

        % Scenario 1: Biting rates differ based on HABAC/MABAC
        % Scenario 0: Biting rates do not differ
        % Each line below assumes the biting rate is the same for 7 days at
        % a time; the biting rate for week 5 is used for week 5 and beyond
        if biting_vary ==1
            % Uninfected Biting rates for HABAC
            p.cc(1:7) = p.c*biting_mat(1,1);
            p.cc(8:14) = p.c*biting_mat(1,2);
            p.cc(15:21) = p.c*biting_mat(1,3);
            p.cc(22:28) = p.c*biting_mat(1,4);
            p.cc(29:p.AGES) = p.c*biting_mat(1,5);

            % Uninfected Biting rates for MABAC
            p.cc2(1:7) = p.c*biting_mat(2,1);
            p.cc2(8:14) = p.c*biting_mat(2,2);
            p.cc2(15:21) = p.c*biting_mat(2,3);
            p.cc2(22:28) = p.c*biting_mat(2,4);
            p.cc2(29:p.AGES) = p.c*biting_mat(2,5);

            % Infected Biting Rates for HABAC
            p.cci(1:7) = p.c*biting_mat_inf(1,1);
            p.cci(8:14) = p.c*biting_mat_inf(1,2);
            p.cci(15:21) = p.c*biting_mat_inf(1,3);
            p.cci(22:28) = p.c*biting_mat_inf(1,4);
            p.cci(29:p.AGES) = p.c*biting_mat_inf(1,5);
            
            % Infected Biting Rates for MABAC
            p.cc2i(1:7) = p.c*biting_mat_inf(2,1);
            p.cc2i(8:14) = p.c*biting_mat_inf(2,2);
            p.cc2i(15:21) = p.c*biting_mat_inf(2,3);
            p.cc2i(22:28) = p.c*biting_mat_inf(2,4);
            p.cc2i(29:p.AGES) = p.c*biting_mat_inf(2,5);
        else
            % Because biting rates do not differ for HABAC/MABAC
            % mosquitoes, we set the rate for uninfected mosquitoes and
            % infected mosquitoes here to be the same for both HABAC and
            % MABAC mosquitoes.
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

        % Set up betaVH rates depending on the scenario - we do not
        % consider this in the present study, but this code allows for
        % exploring this question
        % Scenario 1: HABAC and MABAC rates differ
        % Scenario 0: they do not differ 
        if betaVH_vary==1
            p.betaVH = betaVH0;
            p.betaVH2 = betaVH02;
        else
            p.betaVH = p.betaVH0;
            p.betaVH2 = p.betaVH0;
        end
    
        % Set up betaHV rates depending on the scenario
        % Scenario 1: HABAC and MABAC rates differ based on P. yoelii
        % Secnario 2: HABAC and MABAC rates do not differ, based on P. falciparum
        % Scenario 3: HABAC and MABAC rates differ, based on P. falciparum
        % Scenario 0: HABAC and MABAC rates do not differ, based on P. yoelii 
        if betaHV_vary==1           % P. yoelii with MABAC
            p.betaHV = betaHV0;
            p.betaHV2 = betaHV02;
        elseif betaHV_vary==2       % P. falciparium without MABAC
            p.betaHV = betaHV0f;
            p.betaHV2 = betaHV0f;
        elseif betaHV_vary==3       % P. falciparium with MABAC
            p.betaHV = betaHV0f;
            p.betaHV2 = betaHV02f;
        else                        % P. yoelii without MABAC
            p.betaHV = betaHV0;
            p.betaHV2 = betaHV0;
        end
    
    
        % Fecundity Values
        
        f_choice = 7;       % HABAC
        f_choice2 = 6;      % MABAC


        % Setup the fraction of females who lay eggs and the average number
        % of eggs laid per day per female

        p.frac_fecund = zeros(1,p.AGES);
        p.no_eggs = zeros(1,p.AGES);
        p.frac_fecund2 = p.frac_fecund;
        p.no_eggs2 = p.no_eggs;
        

        % Vectors are Divided by ages 2-6 days, 7-11 days, and older
        
        % HABAC numbers
        p.frac_fecund(2:6) = frac_fecund_mat(f_choice,1);
        p.frac_fecund(7:11) = frac_fecund_mat(f_choice,2);
        p.frac_fecund(12:p.AGES) = frac_fecund_mat(f_choice,3);
        
        p.no_eggs(2:6) = no_eggs_mat(f_choice,1)/2;
        p.no_eggs(7:11) = no_eggs_mat(f_choice,2)/2;
        p.no_eggs(12:p.AGES) = no_eggs_mat(f_choice,3)/2;
    
        % MABAC Nubmers for Scenario 1 in which values differ
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
    
      
    
        % Blood feeding Tendency for uninfected mosquitoes
        
        % Set up Blood feeding tendency vectors
        % Vectors will be organized by age in days
        p.av = zeros(p.AGES,1);
        p.av2 = p.av;
    
        % Second Blood Meal at Day 4
        p.av(1) = 0;
        
        % Different possible distributions for bloodmeals to be spread
        % across 3 days

        % dist = [.25 .5 .25];
        % dist = [0 1 0];
        dist = [0.333 0.333 0.333];
    
        % Index of the blood meal tendency vectors
        bf_num = 7; % HABAC
        bf_num2 = 6; % MABAC
    
        % first 4 bloodmeals
        p.av(2:4) = blood_feeding4(bf_num)*dist;    
        p.av(7:9) = 0.5*(blood_feeding4(bf_num)+blood_feeding14(bf_num))*dist;
        p.av(12:14) = blood_feeding14(bf_num)*dist;
        p.av(17:19) = .778*p.wklyfeed(3)*dist;

        % 5th and beyond
        p.av(22:24) = .778*p.wklyfeed(4)*dist;
        p.av(27:29) = .778*p.wklyfeed(4)*dist;
        p.av(32:34) = .778*p.wklyfeed(4)*dist;
        p.av(37:39) = .778*p.wklyfeed(4)*dist;
        p.av(42:44) = .778*p.wklyfeed(4)*dist;
        p.av(47:49) = .778*p.wklyfeed(4)*dist;

    
        % If uninfected MABAC and HABAC mosquitoes have differenent tendencies,
        % bf_vary=1; otherwise bf_Vary=0
      if bf_vary==1
    
        % first 4 bloodmeals
        p.av2(2:4) = blood_feeding4(bf_num2)*dist;
        p.av2(7:9) = 0.5*(blood_feeding4(bf_num2)+blood_feeding14(bf_num2))*dist;
        p.av2(12:14) = blood_feeding14(bf_num2)*dist;
        p.av2(17:19) = .778*p.wklyfeed(3)*dist;

        % 5th and beyond
        p.av2(22:24) = .778*p.wklyfeed(4)*dist;
        p.av2(27:29) = .778*p.wklyfeed(4)*dist;
        p.av2(32:34) = .778*p.wklyfeed(4)*dist;
        p.av2(37:39) = .778*p.wklyfeed(4)*dist;
        p.av2(42:44) = .778*p.wklyfeed(4)*dist;
        p.av2(47:49) = .778*p.wklyfeed(4)*dist;       
      else 
          p.av2 = p.av;    
      end
    
       % Blood feeding Tendency for uninfected mosquitoes
       bfi_num = 7;     %HABAC
       bfi_num2 = 6;    %MABAC
       
       p.avi = p.av;
       p.avi2 = p.av;
       
       % First 4 bloodmeals
       p.avi(2:4) = blood_feeding4inf(bfi_num)*dist;
       p.avi(7:9) = (0.5)*(blood_feeding4inf(bfi_num)+blood_feeding11inf(bfi_num))*dist;
       p.avi(12:14) = blood_feeding11inf(bfi_num)*dist;
       p.avi(17:19) = .778*p.wklyfeed(3)*dist;

       % 5th bloodmeal and beyond
       p.avi(22:24) = .778*p.wklyfeed(4)*dist;
       p.avi(27:29) = .778*p.wklyfeed(4)*dist;       
       
       % If infected MABAC and HABAC mosquitoes have differenent tendencies,
       % bf_vary_inf=1; otherwise bf_vary_inf=0
       if bf_vary_inf==1
           % First 4 bloodmeals
           p.avi2(2:4) = blood_feeding4inf(bfi_num2)*dist;
           p.avi2(7:9) = (0.5)*(blood_feeding4inf(bfi_num2)+blood_feeding11inf(bfi_num2))*dist;
           p.avi2(12:14) = blood_feeding11inf(bfi_num2)*dist;           
           p.avi2(17:19) = .778*p.wklyfeed(3)*dist;

           % Fifth bloodmeal and beyond
           p.avi2(22:24) = .778*p.wklyfeed(4)*dist;
           p.avi2(27:29) = .778*p.wklyfeed(4)*dist;
       else
            p.avi2 = p.avi;
       end
       
        
    
        %% Storage Setup for model output
        % Type 1 is HABAC (mosquitoes) or non-severe malaria (humans)
        % Type 2 is MABAC (mosquitoes) or severe malaria (humans)
        
        % Human Compartments
        % 1-5 are S, E, I1, I2, R, respectively
        H = zeros(5, p.DAYS);
        
        % Juvenile Vectors
        JVEC = zeros(1, p.DAYS);
        
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
        
        % Total Humans, Total Vectors
        NH = zeros(1,p.DAYS);
        NV = zeros(1,p.DAYS);
        
        % Store Infected vectors, forces of infections
        INFV = zeros(1,p.DAYS);
        lambdaVH = zeros(1,p.DAYS);
        lambdaHV11 = zeros(1,p.DAYS);
        lambdaHV12 = zeros(1,p.DAYS);

        INFV2 = zeros(1,p.DAYS);
        lambdaVH2 = zeros(1,p.DAYS);
        lambdaHV21 = zeros(1,p.DAYS);
        lambdaHV22 = zeros(1,p.DAYS);

        %% Run Model for Burn-in Period
        
        % Initialize Human Population
        NH0 = p.HumanPop;
        H(:,1) = [NH0 0 0 0 0];

        % Initialize Vector Population
        SV(1,1) = 100;
        JV(:,1) = 100;    

        p.DAYS0=p.DAYS;

        % Run burn-in model to allow vector population to equilibrate
        output0 = discrete_eqns(p, NH, NV, H, SV, SV2, EV, IV, EV2, IV2, INFV, INFV2, lambdaHV11, lambdaHV12, lambdaHV21, lambdaHV22, lambdaVH, lambdaVH2, EGGS, JV);

        % Update vector population at equilibrium
        JV(:,1) = output0.JV(:,end);
        SV(:,1) = output0.SV(:,end);    
        
        % Introduce new human infections 
        H(3,1) = p.IntroCase;
        H(1,1) = H(1,1)-H(3,1);
      
        % Display the combination of scenario and sev_mal vector value
        display([jj ii])
        
        % Run Model        
        output = discrete_eqns(p, NH, NV, H, SV, SV2, EV, IV, EV2, IV2, INFV, INFV2, lambdaHV11, lambdaHV12, lambdaHV21, lambdaHV22, lambdaVH, lambdaVH2, EGGS, JV);


        % Store outputs for humans for each simulation - Type 1
        prevalence(jj,:,ii) = output.H(3,:);                % Total cases each day
        incidence(jj,:,ii) = output.new_cases;              % Total new cases each day
        total_cases(jj,:,ii) = cumsum(output.new_cases);    % cumulative cases each day 
        
        % Store outputs for humans for each simulation - Type 2
        prevalence2(jj,:,ii) = output.H(4,:);               
        incidence2(jj,:,ii) = output.new_cases2;
        total_cases2(jj,:,ii) = cumsum(output.new_cases2);
    
      

    end
end



%% Calculate Peak and Timing of Peak of Incidence

tt = 0:(p.DAYS-1);

% Storage vectors
% max_inc: the value of the incidence at its maximum
% time_max: the timing of the max incidence
% duration: the total duration of the outbreak 
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

%% Save Output as a .mat file to use to generate figures 
save(fileName) 


%% Plot for Paper 
% Please note that for visibility, some plot bounds may need to be
% modified depending on which scenarios are being run. 

% This line chooses the color scheme for curve plots. 
xCols = lines;

% a and b set x limit values for panel a
a = 150;    b = 3.6*365;

% Sets the time at which the calculations are made for plots b-e.
ToC = 4*365;

% [r g b] Color Options for Bar Plots
bc1 = [.3 .3 .3];
bc2 = [.4 .4 .9];

% X limits for panels b-e for consistency
xLimz = [.2-.15 1.15];

% Figure Parameters

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
temp = squeeze(total_cases(:,ToC,:)+total_cases2(:,ToC,:));

prctChangeC = (temp(2,:)-temp(1,:))./(temp(1,:))*100;
b1 = bar(ax3,sev_vec,prctChangeC,'facecolor','flat');

ylabel(ax3,{'% change in total cases within';'four years of introduction'})

xlabel(ax3,'fraction of cases of severe malaria (\delta)')

ylim(ax3,[0 70])
xlim(ax3,xLimz)
b1(1).CData = bc1;

set(ax3,'fontsize',14)
ax3; hold off;

xx1 = get(ax3,'XLim');
yy1 = get(ax3,'YLim');x1x = xx1(1)+.05*(xx1(2)-xx1(1));
y1y = yy1(1)+.9*(yy1(2)-yy1(1));

text(ax3,x1x,y1y,'(b)','FontSize',14)

%%%%%%% PANEL (C) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate % Change in size of outbreak peak from HABAC to including MABAC effects
prctChangeP = (max_inc(2,:)-max_inc(1,:))./(max_inc(1,:))*100;

b1 = bar(ax4,sev_vec,prctChangeP,'facecolor','flat');

ylabel(ax4,{'% change in size';'of outbreak peak'})
xlabel(ax4,'fraction of cases of severe malaria (\delta)')
ylim(ax4,[0 120])
xlim(ax4,xLimz)
set(ax4,'fontsize',14)

b1(1).CData = bc1;

xx1 = get(ax4,'XLim');
yy1 = get(ax4,'YLim');

x1x = xx1(1)+.05*(xx1(2)-xx1(1));
y1y = yy1(1)+.9*(yy1(2)-yy1(1));

text(ax4,x1x,y1y,'(c)','FontSize',14)


%%%%%%% PANEL (D) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate % Change in timing of outbreak peak from HABAC to including MABAC effects
prctChangeT = (time_max(2,:)-time_max(1,:))./(time_max(1,:))*100;

b1 = bar(ax5,sev_vec,prctChangeT,'facecolor','flat');
ylabel(ax5,{'% change in time';'of outbreak peak'})
xlabel(ax5,'fraction of cases of severe malaria (\delta)')
set(ax5,'fontsize',14)

ylim(ax5,[-50 10])
xlim(ax5,xLimz)

b1(1).CData = bc1;

xx1 = get(ax5,'XLim');
yy1 = get(ax5,'YLim');

x1x = xx1(1)+.05*(xx1(2)-xx1(1));
y1y = yy1(1)+.9*(yy1(2)-yy1(1));

text(ax5,x1x,y1y,'(d)','FontSize',14)

%%%%%%% PANEL (E) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate % Change in endemic prevalence from HABAC to including MABAC effects

temp = squeeze(prevalence(:,end-1,:)+prevalence2(:,end-1,:));
prctChangeE = (temp(2,:)-temp(1,:))./(temp(1,:))*100;
b1 = bar(ax6,sev_vec,prctChangeE,'facecolor','flat');


ylabel(ax6,{'% change in';'endemic prevalence'})

xlabel(ax6,'fraction of cases of severe malaria (\delta)')
ylim(ax6,[0 30])
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


% saveas(gcf,['Modeling Paper/' fileName '.pdf'])
saveas(gcf,[fileName '.pdf'])

%% Model Equations 
function output = discrete_eqns(p, NH, NV, H, SV, SV2, EV, IV, EV2, IV2, INFV, INFV2, lambdaHV11, lambdaHV12, lambdaHV21, lambdaHV22, lambdaVH, lambdaVH2, EGGS, JV)

    % Storage for new cases
    new_cases = zeros(1,p.DAYS);
    new_cases2 = new_cases;

    for t=1:(p.DAYS-1)
    
        % Calculate Total Human Population each day - this should not
        % change for the present paper
        NH(t) = sum(H(:,t));

        % Calculate the total vector population for each day
        NV(t) = sum(SV(:,t)+EV(:,t)+IV(:,t) + SV2(:,t) + IV2(:,t) + EV2(:,t));
        
        % Calculate total Type1 malaria Cases and multiply times biting rates and
        % tendency to bite 
        INFV(t) = sum(p.cci.*p.avi.*IV(:,t)); 

        % Daily force of infection for Type 1
        lambdaVH(t) = 1-exp(-INFV(t)/NH(t));

        % Calculate total type 2 malaria cases and multiply times biting rates and
        % tendency to bite
        INFV2(t) = sum(p.cc2i.*p.avi2.*IV2(:,t));    

        % Daily force of infection for Type 2
        lambdaVH2(t) = 1-exp(-INFV2(t)/NH(t));
    
        % Equations for Human Dynamics
        H(1,t+1) = H(1,t)*(1-p.betaVH*lambdaVH(t)-p.betaVH2*lambdaVH2(t)) + p.psiH*H(5,t);
        H(2,t+1) = H(2,t)*(1-p.sigmaH) + (p.betaVH*lambdaVH(t)+p.betaVH2*lambdaVH2(t))*H(1,t);
        H(3,t+1) = H(3,t)*(1-p.gammaH) + (1-p.sevMal)*p.sigmaH*H(2,t);
        H(4,t+1) = H(4,t)*(1-p.gammaH2) + (p.sevMal)*p.sigmaH*H(2,t);
        H(5,t+1) = H(5,t)*(1-p.psiH) + p.gammaH*H(3,t) + p.gammaH2*H(4,t);
        
        % Calculate New Cases (non-severe)
        new_cases(t+1) = (1-p.sevMal)*p.sigmaH*H(2,t);
        % Calculate New Severe Cases 
        new_cases2(t+1) = p.sevMal*p.sigmaH*H(2,t);
       
        % Loop to run through ages for mosquitoes     
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
    
        % New Eggs
        EGGS(t) = EGGS(t)+p.frac_fecund(a)*p.no_eggs(a)*(SV(a,t) + EV(a,t) + IV(a,t)) ...
            + p.frac_fecund2(a)*p.no_eggs2(a)*(SV2(a,t) + EV2(a,t) + IV2(a,t));
    
    end
    
    
        % Juvenile and Newly Emerged Vector Dynamics        
        JV(1,t+1) = p.muE*JV(1,t)*(1-p.hatch) + EGGS(t);
        JV(2,t+1) = p.muE*p.hatch*JV(1,t) + p.muJV*exp(-((p.a*JV(2,t))^p.b))*JV(2,t)*(1 - p.nuJV);        
        SV(1, t+1) = (1/2)*p.nuJV*p.muJV*exp(-((p.a*JV(2,t))^p.b))*JV(2,t);
    
    end

    % Output Storage
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