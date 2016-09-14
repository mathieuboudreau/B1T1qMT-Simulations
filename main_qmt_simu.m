%% MTSPGR QA: Stage 2
% Version: 0.3
%
% PROJECT MANAGER: Mathieu Boudreau (_MJB)
% DESCRIPTION: This code investigates the error in fitting parameters that
% results in a range T1 observed value different from the true value.
%
%
% LAST EDITED BY: Mathieu Boudreau (_MJB)
% LAST MODIFIED ON: August 1st 2012
%


%% Clear Matlab
%

clear all
clc
close all

%% Startup
%

%%% Start timing %%%
tic
%%%

startup;                    

outputDir='/Users/mathieuboudreau/Work/Projects/B1T1qMT_simu/Stage2_Opt/';
mkdir([outputDir 'b0'])
mkdir([outputDir 'b1'])
mkdir([outputDir 'qmt'])
mkdir([outputDir 't1'])

disp(outputDir)

%% Set Measurement Parameters
%

%MeasurementParameters = 'experimental_protocols/sled_expt';
%MeasurementParameters = 'CustomProtocol1.5T.mat';
%MeasurementParameters = 'OptimumProtocol3T_10p';
MeasurementParameters = 'UKProtocol_3T';

%protocolFlag = 'sled';
protocolFlag = 'custom';

% Set dimensions, resolution
vol_dims= [151 101];
vol_resolution=[1 1 1];

%% Set Tissue Parameters
%

TissueParameters = 'wm_mt_system_params_3t';
%TissueParameters = 'WM_tissue_parameters_3T';

%% Run simulation
%

load(MeasurementParameters)

switch protocolFlag
   case 'custom'
       [sled_flips sled_offsets sled_pulse_durations sled_trs sled_excite_flips] = create_protocol_sled(protocol);
end
for ii=1:length(protocol)
   mt_measurement(ii) = sim_qmtspgr_singlepoint(protocol,ii,TissueParameters, protocolFlag);
end

mt_measurement=cell2mat(mt_measurement');

%% Produce qMT, B0, B1, T1 and dR1 Minc files
%

%STAGE 0
%b1t1qmtSimu_Stage0(outputDir,mt_measurement,vol_dims,vol_resolution, protocol, sled_trs, sled_excite_flips) % B1 and T1 images are perpendicular 2D gradients for their respective ranges


%STAGE 1
%b1t1qmtSimu_Stage1(outputDir,mt_measurement,vol_dims,vol_resolution, protocol, sled_trs, sled_excite_flips) % B1 and T1 images are perpendicular 2D gradients for their respective ranges

%STAGE 2
%b1t1qmtSimu_Stage2(outputDir,mt_measurement,vol_dims,vol_resolution, protocol, sled_trs, sled_excite_flips) % B1 is a 2D gradient of its ranges, T1 is a perpendicular gradient of VFA calculated T1 (from B1), with a range of T1 values perpendicular to the gradient of B1

%STAGE 2singleslice
b1t1qmtSimu_Stage2_singleline(outputDir,mt_measurement,vol_dims,vol_resolution, protocol, sled_trs, sled_excite_flips) % B1 is a 2D gradient of its ranges, T1 is a perpendicular gradient of VFA calculated T1 (from B1), with a range of T1 values perpendicular to the gradient of B1


%STAGE 3
%b1t1qmtSimu_Stage3(outputDir,mt_measurement,vol_dims,vol_resolution,
%protocol, sled_trs, sled_excite_flips) % Calculates noisy MT measurements, for noise sensitivity analysis of "B1 insensitivity" 


%% ROI
%

roi_vol=zeros(vol_dims);
%roi_vol=zeros(vol_dims(2),vol_dims(1)); % For stage3, I know it's shitty, I should fix that it's the wrong order soon, but didn't want to bother it while working.

roi_vol(51,51)=1;

ROI = minc;
ROI.fileName=[outputDir 'roi.mnc'];
ROI=ROI.setVolume(roi_vol);
ROI=ROI.setResolution(vol_resolution);
ROI.saveMinc;
%% End Program
%

%%% End timing %%%
howLong = toc;
time_result = ['Elapsed time is ', num2str((howLong-rem(howLong,60))/60), ' minutes and ' num2str(rem(howLong,60)) ' seconds'];
disp(time_result)
%%%