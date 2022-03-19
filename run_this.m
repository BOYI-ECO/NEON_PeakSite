
for site = 1:14

    site

    % set work space
    cd('C:\Users\boyi\Desktop\8.MCMC_Uniform_ForLabIcubation_LoopSite_Peak\20200210_dream_zs_restart_with_SOM');

    % Prepare obser_data.csv for each site %
    Obser_data = readmatrix('.\OBSER_data.xlsx');  % Direct read model input %
    %Obser_data = xlsread('.\OBSER_data.xlsx');  % Direct read model input %
    col_name = {' ','temp','water_cont','water_sat','water_fc','Pclay','AOM1_C','AOM1_N','AOM2_C','AOM2_N','AOM3_C','AOM3_N','SMB1_C','SMB1_N','SMB2_C','SMB2_N','SMR_C','SMR_N','NOM_C','NOM_N','MOM_C','MOM_N','NH4','NO3','litter_C','litter_N','litter_lignin','soil_C(mg)','soil_N(mg)','Litter_C(mg)','Litter_N(mg)'};
    %col_name = {'name','temp','water_cont','water_sat','water_fc','Pclay','AOM1_C','AOM1_N','AOM2_C','AOM2_N','AOM3_C','AOM3_N','SMB1_C','SMB1_N','SMB2_C','SMB2_N','SMR_C','SMR_N','NOM_C','NOM_N','MOM_C','MOM_N','NH4','NO3','litter_C','litter_N','litter_lignin','soil_C','soil_N','Litter_C','Litter_N'};
    obser_data = Obser_data(site,:);
    obser_data(1) = 1;
    obser_data = array2table(obser_data);
    obser_data.Properties.VariableNames(:) = col_name;
    writetable(obser_data,'.\Rh_1\obser_data.csv')

    % Get MCMC setting  
    N = 12;                  % Number of parallel chains
    T = 600;                 % Number of iterations
    Npar = 33;               % Dimension of the model parameters
    % Nout = 350;            % Dimension of the model responses
    copyexample(N);          % parallel 
    
    % Get prior info
    xmin = importdata('para_min.txt')';      % Lower bound of the prior model parameters
    xmax = importdata('para_max.txt')';      % Upper bound of the prior model parameters
    para_idx = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33];
    xmin = xmin(para_idx);
    xmax = xmax(para_idx);
    range = [xmin' xmax'];                   % Range of the prior model parameters
    tic;

    % Method1: Lab Obervation data %
    lab_obser = readmatrix('C:\Users\boyi\Desktop\8.MCMC_Uniform_ForLabIcubation_LoopSite_Peak\20200210_dream_zs_restart_with_SOM\lab_obser.csv');  % Direct read Observation Data %
    %lab_obser = xlsread('.\lab_obser.csv');  % Direct read Observation Data %
    yreal     = lab_obser(:,site);            % Second Column is CO2 emission data % 
    Nout      = size(yreal,1);                % Total numbers of Output comparision %
    sd        = yreal*0.1;
    Obs       = yreal + sd.*randn(size(yreal));

    % Method2: Data Comparison %
    % xreal = prior(1,Npar,xmin,xmax);        % The reference parameters drawn from the prior distribution
    xreal = importdata('true_para.txt');       % Just used for drawing the prior(true) para in output graph 
    xreal = xreal(para_idx);                  
    % yreal = forwardmodel(xreal,1);
    % sd = yreal*0.1;                         % Standard deviation of the measurement errors
    % Obs = yreal + sd.*randn(size(yreal));   % The measurements

    % Prior Space %
    m0 = max(N,60*Npar);  %% (T*N > *Npar)When N=8,T=500, and use 200*Npar: Error occur (Dont hnow why)!
    Z0 = prior(m0,Npar,xmin,xmax);
    Z_prior = Z0;
    
    for i = 1:5 
        i
        [x,fx] = dream_zs(N,T,Npar,Nout,Obs,sd,range,Z0);
        chain  = GenParSet(x);
        % drawplot(chain,xreal,range,2,3);
        Z0 = chain(end-m0+1:end,:);    % if error occur here, T/N/X(*Npar) is not correct  
        % Z0 = chain(randsample(size(chain,1),m0),:);
    end
    
    outputfolder = 'C:\Users\boyi\Desktop\8.MCMC_Uniform_ForLabIcubation_LoopSite_Peak\20200210_dream_zs_restart_with_SOM';
    %Site_name = {'BONA1','BONA2','BONA3','BONA4','CPER1','CPER2','CPER3','CPER4','DSNY1','DSNY2','DSNY3','DSNY4','GRSM1','GRSM2','GRSM3','GRSM4','HARV1','HARV2','HARV3','HARV4','KONZ1','KONZ2','KONZ3','KONZ4','LENO1','LENO2','LENO3','LENO4','NIWO1','NIWO2','NIWO3','NIWO4','ONAQ1','ONAQ2','ONAQ3','ONAQ4','OSBS1','OSBS2','OSBS3','OSBS4','PUUM1','PUUM2','PUUM3','PUUM4','SCBI1','SCBI2','SCBI3','SCBI4','SJER1','SJER2','SJER3','SJER4','SRER1','SRER2','SRER3','SRER4','TALL1','TALL2','TALL3','TALL4','TOOL1','TOOL2','TOOL3','TOOL4','UNDE1','UNDE2','UNDE3','UNDE4','WOOD1','WOOD2','WOOD3','WOOD4','WREF1','WREF2','WREF3','WREF4','YELL1','YELL2','YELL3','YELL4'};
    Site_name = {'GRSM3','GRSM4','HARV1','HARV4','LENO1','LENO2','LENO3','LENO4','OSBS1','OSBS2','SCBI1','UNDE3','WREF2','WREF4'};

    file_name = [Site_name{site},'.mat'];
    save (file_name)

    toc
    copyexample(N,-1);

    %% plot posterior density of paras
    %for i = 1:Npar
    %    subplot(5,6,i)
    %    [k,value] = ksdensity(Z0(:,i));
    %    plot(value,k,'Linewidth',2)
    %    hold on;
    %    [k_p,value_p] = ksdensity(Z_prior(:,i));
    %    plot(value_p,k_p,'Linewidth',2,'color','k');
    %    hold on;
    %    xreal = importdata('true_para.txt');
    %    xreal = xreal(para_idx);
    %    plot([xreal(i) xreal(i)],[0 max([k,k_p])*1.1],'LineWidth',2,'Color',[1 0 0]);
    %    ylim([0,max([k,k_p])*1.1]);
    %end
    
    %file_name = [Site_name{site},'.png'];
    %saveas(gcf,file_name);

    % clear all
    clear
    tic

end

% 1) cal_Rh.py file: "day_selected" and "column[;,-3:]" for comparision


