function [Stage1,Stage2,Stage3] = analysis__systemMatrix(dataDir,vol_dims_allStages)
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
Stage1.R1f = reshape(fit.R1(:,1),Stage1.vol_dims);
Stage1.R1r = reshape(fit.R1(:,2),Stage1.vol_dims);

Stage1.T2f = reshape(fit.T2(:,1),Stage1.vol_dims);
Stage1.T2r = reshape(fit.T2(:,2),Stage1.vol_dims);

%Load Errors data
Stage1.dF = reshape(fit.df,Stage1.vol_dims);
Stage1.E2 = reshape(fit.e2,Stage1.vol_dims);
Stage1.dF(Stage1.dF>1)=0;
Stage1.E2(Stage1.E2>10^-3)=0;

Stage1.F_Perc=((Stage1.F-trueF)./trueF.*100);
Stage1.kf_Perc=((Stage1.kf-truekf)./truekf.*100);

Stage1.B1_range=Stage1.B1(1,:);
Stage1.T1_range=Stage1.T1(:,1)';

for ii=1:Stage1.vol_dims(1)
    for jj=1:Stage1.vol_dims(2)
        system_matrix_rf_off{ii,jj}=[-(Stage1.R1f(ii,jj)+Stage1.kf(ii,jj)), Stage1.kf(ii,jj)./Stage1.F(ii,jj) ; Stage1.kf(ii,jj), -Stage1.kf(ii,jj)./Stage1.F(ii,jj) ];
        eig_system_matrix_rf_off{ii,jj}=eig([-(Stage1.R1f(ii,jj)+Stage1.kf(ii,jj)), Stage1.kf(ii,jj)./Stage1.F(ii,jj) ; Stage1.kf(ii,jj), -Stage1.kf(ii,jj)./Stage1.F(ii,jj) ]);
        eig_system_matrix_rf_off_EIG_1(ii,jj)=eig_system_matrix_rf_off{ii,jj}(1);
        eig_system_matrix_rf_off_EIG_2(ii,jj)=eig_system_matrix_rf_off{ii,jj}(2);
        
    end
end
load('/Users/mathieuboudreau/Work/mboudrea-phd-code/MATLAB/MyCollections/b1t1qmt_simu/parameters/pulse_sequence/qmt/UKProtocol_3T.mat')

%gamma= 267.513*10^6;
theta_nom_1=protocol(1,1)*pi/180;
pulse_duration=protocol(1,3);

offres_freq=protocol(1,2);
omega_rms = theta_nom_1./pulse_duration;

R1r=1;

dR1f_dkf = 1-( (R1r-1./trueT1)./ ((R1r-1./trueT1)+truekf./trueF)  ) - truekf./trueF.* ( (R1r-1./trueT1)./ ((R1r-1./trueT1)+truekf./trueF)  ) .* ((R1r-1./trueT1)+truekf./trueF).^-1
dR1f_dF = 1- (truekf./trueF).^2.* ( (R1r-1./trueT1)./ ((R1r-1./trueT1)+truekf./trueF)  ) .* ((R1r-1./trueT1)+truekf./trueF).^-1
dR1f_dB1=(Stage1.R1f(80,66)-Stage1.R1f(80,67))./(Stage1.B1(80,66)-Stage1.B1(80,67))


%% For B1(26:76) ie. B1=0.75 to 1.25, polyfit of order 2 gives the following relationship for T1
%  T1_fitted = 2.9962*B1^2 - 7.9601 *B1 +5.8621 seconds

%% Go back to initial directory
%

cd(oldDir)

end

