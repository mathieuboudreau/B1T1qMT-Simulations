function [Stage1,Stage2,Stage3] = analysis__full_project(dataDir,vol_dims_allStages)
%ANALYSIS__FULL_PROJECT Summary of this function goes here
%   Detailed explanation goes here
plot_settings

oldDir=cd;

%% Load Stage 0 data
% This is the "true" fit of the tissue

%Change folder
cd([dataDir '/Stage0'])

% Open data files
load('preparedData.mat')
load('round2.mat')

%Load B1, T1, F data
trueB1 = data.b1;
trueT1 = 1./data.R1obs;
trueF = fit.f;
truekf = fit.kf;

%% Load Stage 1 data
%

%Change folder
cd([dataDir '/Stage1'])

% Open data files
load('preparedData.mat')
load('round2.mat')

%Get Volume Dimensions
Stage1.vol_dims=vol_dims_allStages.Stage1;

%Load B1, T1, F data
Stage1.B1 = reshape(data.b1,Stage1.vol_dims);
Stage1.T1 = reshape(1./data.R1obs,Stage1.vol_dims);
Stage1.F = reshape(fit.f,Stage1.vol_dims);
Stage1.kf = reshape(fit.kf,Stage1.vol_dims);

%Load Errors data
Stage1.dF = reshape(fit.df,Stage1.vol_dims);
Stage1.E2 = reshape(fit.e2,Stage1.vol_dims);
Stage1.dF(Stage1.dF>1)=0;
Stage1.E2(Stage1.E2>10^-3)=0;

Stage1.F_Perc=((Stage1.F-trueF)./trueF.*100);
Stage1.kf_Perc=((Stage1.kf-truekf)./truekf.*100);

Stage1.B1_range=Stage1.B1(1,:);
Stage1.T1_range=Stage1.T1(:,1)';


% Plot preliminary maps
figure(101)
    imagesc(Stage1.B1_range,Stage1.T1_range,Stage1.B1)
    my_xlabel('B_{1} (n.u.)')
    my_ylabel('T_{1} (s)')
    title('Map of B_{1} values used in Stage 1')


figure(102)
    imagesc(Stage1.B1_range,Stage1.T1_range,Stage1.T1)
    my_xlabel('B_{1} (n.u.)')
    my_ylabel('T_{1} (s)')
    title('Map of T_{1} values used in Stage 1')


figure(103)
    imagesc(Stage1.B1_range,Stage1.T1_range,Stage1.F)
    caxis([0,2*trueF])
    my_xlabel('B_{1} (n.u.)')
    my_ylabel('T_{1} (s)')
    title('Map of fitted F values in Stage 1')


figure(104)
    imagesc(Stage1.B1_range,Stage1.T1_range,Stage1.kf)
    caxis([0,2*truekf])
    my_xlabel('B_{1} (n.u.)')
    my_ylabel('T_{1} (s)')
    title('Map of fitted kf values in Stage 1')

% Plot preliminary maps
figure(113)
    imagesc(Stage1.B1_range,Stage1.T1_range,Stage1.dF)
    my_xlabel('B_{1} (n.u.)')
    my_ylabel('T_{1} (s)')
    title('Map of B_{1} values used in Stage 1')

% Plot preliminary maps
figure(120)
    imagesc(Stage1.B1_range,Stage1.T1_range,Stage1.E2)
    my_xlabel('B_{1} (n.u.)')
    my_ylabel('T_{1} (s)')
    title('Map of B_{1} values used in Stage 1')    
    
%% Load Stage 2 data
%
%Change folder
cd([dataDir '/Stage2'])

% Open data files
load('preparedData.mat')
load('round2.mat')

%Get Volume Dimensions
Stage2.vol_dims=vol_dims_allStages.Stage2;

%Load B1, T1, F data
Stage2.B1 = reshape(data.b1,Stage2.vol_dims);
Stage2.T1 = reshape(1./data.R1obs,Stage2.vol_dims);
Stage2.F = reshape(fit.f,Stage2.vol_dims);
Stage2.kf = reshape(fit.kf,Stage2.vol_dims);

Stage2.F_Perc=((Stage2.F-trueF)./trueF.*100);
Stage2.kf_Perc=((Stage2.kf-truekf)./truekf.*100);


% *** Fix the hardcoded 51 value...? MJB
Stage2.B1_range=Stage2.B1(51,:);
Stage2.T1_range=Stage2.T1(:,51)';

% Plot preliminary maps
figure(201)
    imagesc(Stage2.B1_range,Stage2.T1_range,Stage2.B1)
    my_xlabel('B_{1} (n.u.)')
    my_ylabel('T_{1} (s)')
    title('Map of B_{1} values used in Stage 2')


figure(202)
    imagesc(Stage2.B1_range,Stage2.T1_range,Stage2.T1)
    my_xlabel('B_{1} (n.u.)')
    my_ylabel('T_{1} (s)')
    title('Map of VFA fitted T_{1} values used in Stage 2')


figure(203)
    imagesc(Stage2.B1_range,Stage2.T1_range,Stage2.F), caxis([0,2*trueF])
    my_xlabel('B_{1} (n.u.)')
    my_ylabel('T_{1} (s)')
    title('Map of fitted F values used, using VFA fitted T_{1} values, in Stage 2')


figure(204)
    imagesc(Stage2.B1_range,Stage2.T1_range,Stage2.kf), caxis([0,2*truekf])
    my_xlabel('B_{1} (n.u.)')
    my_ylabel('T_{1} (s)')
    title('Map of fitted kf values used, using VFA fitted T_{1} values, in Stage 2')


%% Load Stage 3 data
%

%Change folder
cd([dataDir '/Stage3'])


for numSteps=1:10
    % Open data files
    load(['round2_' num2str(numSteps) '.mat'])
    if numSteps==1
        Stage3.F=reshape(fit.f,100,101);
        Stage3.kf=reshape(fit.kf,100,101);

    else
        Stage3.F=cat(1,Stage3.F,reshape(fit.f,100,101));
        Stage3.kf=cat(1,Stage3.kf,reshape(fit.kf,100,101));
    end
end

% Open data files
load('preparedData_1.mat')
Stage3.B1=reshape(data.b1,100,101);

Stage3.meanF=mean(Stage3.F);
Stage3.stdF=std(Stage3.F);

Stage3.meankf=mean(Stage3.kf);
Stage3.stdkf=std(Stage3.kf);


figure(301)
    plot(Stage3.B1(51,:),(Stage3.meanF),'LineWidth',4)
    hold on
    plot(Stage3.B1(51,:),Stage2.F(51,:),'r','LineWidth',4)
    errorbar(Stage3.B1(51,1:10:end),Stage3.meanF(1:10:end),Stage3.stdF(1:10:end),'.','LineWidth',4)
    axis([0.49 1.51 0 0.3])
    title('Monte carlo (n=1000) qMT fitting vs Noisless qMT fitting of the pool size ratio F, for a range of B_{1} error')
    legendCell{1}='Noisy qMT - Mean F';
    legendCell{2}='Noiseless qMT - F';
    my_legend(legendCell)

figure(302)
    plot(Stage3.B1(51,:),(Stage3.stdF./Stage3.meanF),'LineWidth',4) % Coefficient of Variation
    title('Coefficient of variation of Monte carlo (n=1000) qMT fitting of F ')
    
figure(303)
    plot(Stage3.B1(51,:),(Stage3.meanF-trueF)./trueF*100,'LineWidth',4)
    hold on
    plot(Stage3.B1(51,:),(Stage2.F(51,:)-trueF)./trueF*100,'r','LineWidth',4)
    title('Relative error of fitted F vs true F for a range of B_{1} error')
    legendCell{1}='Noisy F';
    legendCell{2}='Noisless F';
    my_legend(legendCell)
    
figure(304)
    plot(Stage3.B1(51,:),(Stage3.meanF-Stage2.F(51,:))./Stage2.F(51,:)*100,'LineWidth',4)
    title('Relative error of mean Noisy F vs Noisless F for a range of B_{1} error')

figure(305)
    [y,x]=hist(Stage3.F(:,51),50);
    plot(x,y,'LineWidth',4)
    title('Histogram of Monte Carlo(n=1000) fitted F for B_{1}=1')
    
    
figure(311)
    plot(Stage3.B1(51,:),(Stage3.meankf),'LineWidth',4)
    hold on
    plot(Stage3.B1(51,:),Stage2.kf(51,:),'r','LineWidth',4)%,errorbar(Stage3.B1(51,1:10:end),Stage3.meankf(1:10:end),Stage3.stdkf(1:10:end),'.','LineWidth',4)
    %axis([0.49 1.51 0 0.3])
    title('Monte carlo (n=1000) qMT fitting vs Noisless qMT fitting of the exchange rate kf, for a range of B_{1} error')
    legendCell{1}='Noisy qMT - Mean kf';
    legendCell{2}='Noiseless qMT - kf';
    my_legend(legendCell)
    
figure(312)
    plot(Stage3.B1(51,:),(Stage3.stdkf./Stage3.meankf),'LineWidth',4) % Coefficient of Variation
    title('Coefficient of variation of Monte carlo (n=1000) qMT fitting of kf ')

figure(313)
    plot(Stage3.B1(51,:),(Stage3.meankf-truekf)./trueF*100,'LineWidth',4)
    hold on
    plot(Stage3.B1(51,:),(Stage2.kf(51,:)-truekf)./truekf*100,'r','LineWidth',4)
    title('Relative error of fitted kf vs true kf for a range of B_{1} error')
    legendCell{1}='Noisy kf';
    legendCell{2}='Noisless kf';
    my_legend(legendCell)

figure(314)
    plot(Stage3.B1(51,:),(Stage3.meankf-Stage2.kf(51,:))./Stage2.kf(51,:)*100,'LineWidth',4)
    title('Relative error of mean Noisy kf vs Noisless kf for a range of B_{1} error')

figure(315)
    histKF=Stage3.kf(:,51);
    histKF(histKF>20)=[];
    [y,x]=hist(histKF,50);
    plot(x,y,'LineWidth',4)
    title(['Histogram of Monte Carlo(n=' num2str(length(histKF)) ') fitted F for B_{1}=1. Values above kf=20 were dumped(n=' num2str(1000-length(histKF)) ')'])

%% Miscellanious plots

figure(001)
    imagesc(Stage1.B1_range,Stage1.T1_range,Stage1.F_Perc), caxis([-100,100])
    hold on
    plot(Stage2.B1(51,:), Stage2.T1(51,:),'--k','LineWidth',3)

figure(002)
    imagesc(Stage1.B1_range,Stage1.T1_range,Stage1.F_Perc), caxis([-10,10])
    hold on
    plot(Stage2.B1(51,:), Stage2.T1(51,:),'--k','LineWidth',3)


figure(003)
    plot(Stage2.B1(51,:), Stage2.F_Perc(51,:),'LineWidth',4)
    hold on
    plot(Stage1.B1_range, Stage1.F_Perc(80,:),'r','LineWidth',4)
    
    
figure(004)
    imagesc(Stage1.B1_range,Stage1.T1_range,Stage1.kf_Perc), caxis([-100,100])
    hold on
    plot(Stage2.B1(51,:), Stage2.T1(51,:),'--k','LineWidth',3)

figure(005)
    imagesc(Stage1.B1_range,Stage1.T1_range,Stage1.kf_Perc), caxis([-10,10])
    hold on
    plot(Stage2.B1(51,:), Stage2.T1(51,:),'--k','LineWidth',3)


figure(006)
    plot(Stage2.B1(51,:), Stage2.kf_Perc(51,:),'LineWidth',4)
    hold on
    plot(Stage1.B1_range, Stage1.kf_Perc(80,:),'r','LineWidth',4)
    %axis([0.5 1.5 -20 40])

figure(007)
    plot(Stage2.B1(51,:), Stage2.kf(51,:),'LineWidth',4)
    hold on
    plot(Stage2.B1(51,:), Stage2.kf(41,:),'r','LineWidth',4)
    plot(Stage2.B1(51,:), Stage2.kf(61,:),'k','LineWidth',4)
    
figure(011)
    imagesc(Stage1.B1_range,Stage1.T1_range,Stage1.F_Perc), caxis([-100,100])
    hold on
    plot(Stage2.B1(51,:), Stage2.T1(51,:),'--r','LineWidth',3)
    plot(Stage2.B1(51,:), linspace(0.9,0.9,length(Stage2.B1(51,:))),'--b','LineWidth',3)    

tmp=Stage1.F_Perc(:,67);
tmp(tmp>1000)=1000;
figure(013)
    plot(Stage1.T1_range, tmp,'r','LineWidth',4)
%% Go back to initial directory
%

cd(oldDir)

end

