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


set(0, 'DefaultAxesBox', 'on', 'DefaultAxesLineWidth', 6);
set(0, 'DefaultAxesFontSize', 30, 'DefaultAxesFontWeight', 'bold');


%% Set Measurement Parameters
%

%MeasurementParameters = 'experimental_protocols/sled_expt';
%MeasurementParameters = 'CustomProtocol1.5T.mat';
%MeasurementParameters = 'OptimumProtocol3T_10p';
MeasurementParameters = 'UKProtocol_3T';

%protocolFlag = 'sled';
protocolFlag = 'custom';


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

mt_measurement_true=cell2mat(mt_measurement');

%%
%

load('UKProtocol_3T_b1_5perc_high')

switch protocolFlag
   case 'custom'
       [sled_flips sled_offsets sled_pulse_durations sled_trs sled_excite_flips] = create_protocol_sled(protocol);
end
for ii=1:length(protocol)
   mt_measurement(ii) = sim_qmtspgr_singlepoint(protocol,ii,TissueParameters, protocolFlag);
end

mt_measurement_b1_5perc_high=cell2mat(mt_measurement');

%%
%

MeasurementParameters = 'UKProtocol_3T';
load(MeasurementParameters)
TissueParameters = 'wm_mt_system_params_3t_b1_5perc_high';

switch protocolFlag
   case 'custom'
       [sled_flips sled_offsets sled_pulse_durations sled_trs sled_excite_flips] = create_protocol_sled(protocol);
end
for ii=1:length(protocol)
   mt_measurement(ii) = sim_qmtspgr_singlepoint(protocol,ii,TissueParameters, protocolFlag);
end

mt_measurement_t1_b1_5perc_high=cell2mat(mt_measurement');

%%
%



%%
%

load('UKProtocol_3T_b1_5perc_high')

TissueParameters = 'wm_mt_system_params_3t_b1_5perc_high';

switch protocolFlag
   case 'custom'
       [sled_flips sled_offsets sled_pulse_durations sled_trs sled_excite_flips] = create_protocol_sled(protocol);
end
for ii=1:length(protocol)
   mt_measurement(ii) = sim_qmtspgr_singlepoint(protocol,ii,TissueParameters, protocolFlag);
end

mt_measurement_both=cell2mat(mt_measurement');
%%
%
figure()
semilogx(protocol(1:2:end,2),mt_measurement_true(1:2:end),'s', 'LineWidth',3)
hold on
semilogx(protocol(2:2:end,2),mt_measurement_true(2:2:end),'o', 'LineWidth',3)
my_xlabel('Offset Frequency (Hz)')
my_ylabel('MT Signal (n.u.)')
axis([100 100000 0 1.1])
title('Z-Spectrum for the reference tissue parameters')
legend_cell{1}='MT_FA = 142^{o}';
legend_cell{2}='MT_FA = 426^{o}';
my_legend(legend_cell)

%%
%
figure()
semilogx(protocol(1:2:end,2),(mt_measurement_b1_5perc_high(1:2:end)-mt_measurement_true(1:2:end))./mt_measurement_true(1:2:end).*100,'s', 'LineWidth',3)
hold on
semilogx(protocol(2:2:end,2),(mt_measurement_b1_5perc_high(2:2:end)-mt_measurement_true(2:2:end))./mt_measurement_true(2:2:end).*100,'o', 'LineWidth',3)

semilogx(protocol(1:2:end,2),(mt_measurement_t1_b1_5perc_high(1:2:end)-mt_measurement_true(1:2:end))./mt_measurement_true(1:2:end).*100,'rs', 'LineWidth',3)
semilogx(protocol(2:2:end,2),(mt_measurement_t1_b1_5perc_high(2:2:end)-mt_measurement_true(2:2:end))./mt_measurement_true(2:2:end).*100,'ro', 'LineWidth',3)

my_xlabel('Offset Frequency (Hz)')
my_ylabel('% Difference vs. Reference')
axis([100 100000 -6 6])
title('Relative difference (%) in Z-Spectrums')
legend_cell{1}='Increase in B_1 (142^{o})';
legend_cell{2}='Increase in B_1 (426^{o})';
legend_cell{3}='Decrease in VFA T_1 (142^{o})';
legend_cell{4}='Decrease in VFA T_1 (426^{o})';
my_legend(legend_cell)

%%
%
%semilogx(protocol(:,2),(mean_mt-mt_measurement_true)./mt_measurement_true.*100,'ko')
figure()
semilogx(protocol(1:2:end,2),(mt_measurement_both(1:2:end)-mt_measurement_true(1:2:end))./mt_measurement_true(1:2:end).*100,'s', 'LineWidth',3)
hold on
semilogx(protocol(2:2:end,2),(mt_measurement_both(2:2:end)-mt_measurement_true(2:2:end))./mt_measurement_true(2:2:end).*100,'o', 'LineWidth',3)

my_xlabel('Offset Frequency (Hz)')
my_ylabel('% Difference vs. Reference')
axis([100 100000 -1 4])
title('Relative difference (%) in Z-Spectrums')
legend_cell{1}='Increase in B_1/Decrease in T_1 (142^{o})';
legend_cell{2}='Increase in B_1/Decrease in T_1 (426^{o})';
my_legend(legend_cell)

%%
%
MeasurementParameters = 'UKProtocol_3T';
load(MeasurementParameters)
TissueParameters = 'wm_mt_system_params_3t_F10perc_low';

switch protocolFlag
   case 'custom'
       [sled_flips sled_offsets sled_pulse_durations sled_trs sled_excite_flips] = create_protocol_sled(protocol);
end
for ii=1:length(protocol)
   mt_measurement(ii) = sim_qmtspgr_singlepoint(protocol,ii,TissueParameters, protocolFlag);
end

mt_measurement_F10perc_low=cell2mat(mt_measurement');

%%
%

MeasurementParameters = 'UKProtocol_3T';
load(MeasurementParameters)
TissueParameters = 'wm_mt_system_params_3t_kf10perc_low';

switch protocolFlag
   case 'custom'
       [sled_flips sled_offsets sled_pulse_durations sled_trs sled_excite_flips] = create_protocol_sled(protocol);
end
for ii=1:length(protocol)
   mt_measurement(ii) = sim_qmtspgr_singlepoint(protocol,ii,TissueParameters, protocolFlag);
end

mt_measurement_kf10perc_low=cell2mat(mt_measurement');


%%
%


%%
%

figure()
semilogx(protocol(1:2:end,2),(mt_measurement_F10perc_low(1:2:end)-mt_measurement_true(1:2:end))./mt_measurement_true(1:2:end).*100,'s')
hold on
semilogx(protocol(2:2:end,2),(mt_measurement_F10perc_low(2:2:end)-mt_measurement_true(2:2:end))./mt_measurement_true(2:2:end).*100,'o')



semilogx(protocol(1:2:end,2),(mt_measurement_kf10perc_low(1:2:end)-mt_measurement_true(1:2:end))./mt_measurement_true(1:2:end).*100,'rs')
semilogx(protocol(2:2:end,2),(mt_measurement_kf10perc_low(2:2:end)-mt_measurement_true(2:2:end))./mt_measurement_true(2:2:end).*100,'ro')


%%
%
figure()
semilogx(protocol(1:2:end,2),(mt_measurement_both(1:2:end)-mt_measurement_true(1:2:end))./mt_measurement_true(1:2:end).*100,'s','LineWidth',3)
hold on
semilogx(protocol(2:2:end,2),(mt_measurement_both(2:2:end)-mt_measurement_true(2:2:end))./mt_measurement_true(2:2:end).*100,'o','LineWidth',3)


semilogx(protocol(1:2:end,2),(mt_measurement_kf10perc_low(1:2:end)-mt_measurement_true(1:2:end))./mt_measurement_true(1:2:end).*100,'rs','LineWidth',3)
semilogx(protocol(2:2:end,2),(mt_measurement_kf10perc_low(2:2:end)-mt_measurement_true(2:2:end))./mt_measurement_true(2:2:end).*100,'ro','LineWidth',3)

my_xlabel('Offset Frequency (Hz)')
my_ylabel('% Difference vs. Reference')
axis([100 100000 -1 4])
title('Relative difference (%) in Z-Spectrums')
legend_cell{1}='Increase in B_1/Decrease in T_1 (142^{o})';
legend_cell{2}='Increase in B_1/Decrease in T_1 (426^{o})';
legend_cell{3}='Increase in kf (142^{o})';
legend_cell{4}='Increase in kf (426^{o})';
my_legend(legend_cell)

%%
%

figure()
semilogx(protocol(1:2:end,2),(mt_measurement_both(1:2:end)-mt_measurement_true(1:2:end))./mt_measurement_true(1:2:end).*100,'s','LineWidth',3)
hold on
semilogx(protocol(2:2:end,2),(mt_measurement_both(2:2:end)-mt_measurement_true(2:2:end))./mt_measurement_true(2:2:end).*100,'o','LineWidth',3)


semilogx(protocol(1:2:end,2),(mt_measurement_F10perc_low(1:2:end)-mt_measurement_true(1:2:end))./mt_measurement_true(1:2:end).*100,'rs','LineWidth',3)
semilogx(protocol(2:2:end,2),(mt_measurement_F10perc_low(2:2:end)-mt_measurement_true(2:2:end))./mt_measurement_true(2:2:end).*100,'ro','LineWidth',3)

my_xlabel('Offset Frequency (Hz)')
my_ylabel('% Difference vs. Reference')
axis([100 100000 -1 4])
title('Relative difference (%) in Z-Spectrums')
legend_cell{1}='Increase in B_1/Decrease in T_1 (142^{o})';
legend_cell{2}='Increase in B_1/Decrease in T_1 (426^{o})';
legend_cell{3}='Increase in F (142^{o})';
legend_cell{4}='Increase in F (426^{o})';
my_legend(legend_cell)

%%
%

%%
%
figure()
semilogx(protocol(1:2:end,2),(mt_measurement_b1_5perc_high(1:2:end)-mt_measurement_true(1:2:end))./mt_measurement_true(1:2:end).*100,'s','LineWidth',3)
hold on
semilogx(protocol(2:2:end,2),(mt_measurement_b1_5perc_high(2:2:end)-mt_measurement_true(2:2:end))./mt_measurement_true(2:2:end).*100,'o','LineWidth',3)


semilogx(protocol(1:2:end,2),-(mt_measurement_F10perc_low(1:2:end)-mt_measurement_true(1:2:end))./mt_measurement_true(1:2:end).*100,'rs','LineWidth',3)
semilogx(protocol(2:2:end,2),-(mt_measurement_F10perc_low(2:2:end)-mt_measurement_true(2:2:end))./mt_measurement_true(2:2:end).*100,'ro','LineWidth',3)

my_xlabel('Offset Frequency (Hz)')
my_ylabel('% Difference vs. Reference')
axis([100 100000 -5 1])
title('Relative difference (%) in Z-Spectrums')
legend_cell{1}='Increase in B_1 (142^{o})';
legend_cell{2}='Increase in B_1 (426^{o})';
legend_cell{3}='Decrease in F (142^{o})';
legend_cell{4}='Decrease in F (426^{o})';
my_legend(legend_cell)

%%
%
figure()
semilogx(protocol(1:2:end,2),(mt_measurement_b1_5perc_high(1:2:end)-mt_measurement_true(1:2:end))./mt_measurement_true(1:2:end).*100+(mt_measurement_F10perc_low(1:2:end)-mt_measurement_true(1:2:end))./mt_measurement_true(1:2:end).*100,'s','LineWidth',3)
hold on
semilogx(protocol(2:2:end,2),(mt_measurement_b1_5perc_high(2:2:end)-mt_measurement_true(2:2:end))./mt_measurement_true(2:2:end).*100+(mt_measurement_F10perc_low(2:2:end)-mt_measurement_true(2:2:end))./mt_measurement_true(2:2:end).*100,'o','LineWidth',3)

my_xlabel('Offset Frequency (Hz)')
my_ylabel('% Difference vs. Reference')
axis([100 100000 -5 1])
title('Difference between increase in B_1 and decrease in F')
legend_cell{1}='142^{o}';
legend_cell{2}='426^{o}';
my_legend(legend_cell)


%%
%
figure()
semilogx(protocol(1:2:end,2),(mt_measurement_t1_b1_5perc_high(1:2:end)-mt_measurement_true(1:2:end))./mt_measurement_true(1:2:end).*100+(mt_measurement_F10perc_low(1:2:end)-mt_measurement_true(1:2:end))./mt_measurement_true(1:2:end).*100,'s','LineWidth',3)
hold on
semilogx(protocol(2:2:end,2),(mt_measurement_t1_b1_5perc_high(2:2:end)-mt_measurement_true(2:2:end))./mt_measurement_true(2:2:end).*100+(mt_measurement_F10perc_low(2:2:end)-mt_measurement_true(2:2:end))./mt_measurement_true(2:2:end).*100,'o','LineWidth',3)

my_xlabel('Offset Frequency (Hz)')
my_ylabel('% Difference vs. Reference')
axis([100 100000 0 10])
title('Difference between increase in B_1 and decrease in F')
legend_cell{1}='142^{o}';
legend_cell{2}='426^{o}';
my_legend(legend_cell)


%%
%
%
figure()
semilogx(protocol(1:2:end,2),abs(mt_measurement_b1_5perc_high(1:2:end)-mt_measurement_true(1:2:end))./mt_measurement_true(1:2:end).*100,'s','LineWidth',3)
hold on
semilogx(protocol(2:2:end,2),abs(mt_measurement_b1_5perc_high(2:2:end)-mt_measurement_true(2:2:end))./mt_measurement_true(2:2:end).*100,'rx','LineWidth',3)


my_xlabel('Offset Frequency (Hz)')
my_ylabel('% Difference vs. Reference')
axis([100 100000 -1 4])
title('Increase in B1')
legend_cell{1}='MT FA = 142^{o}';
legend_cell{2}='MT FA = 426^{o}';
my_legend(legend_cell)

figure()
semilogx(protocol(1:2:end,2),abs(mt_measurement_both(1:2:end)-mt_measurement_true(1:2:end))./mt_measurement_true(1:2:end).*100,'s','LineWidth',3)
hold on
semilogx(protocol(2:2:end,2),abs(mt_measurement_both(2:2:end)-mt_measurement_true(2:2:end))./mt_measurement_true(2:2:end).*100,'rx','LineWidth',3)

my_xlabel('Offset Frequency (Hz)')
my_ylabel('% Difference vs. Reference')
axis([100 100000 -1 4])
title('Increase in B1, and change in VFA T1')
legend_cell{1}='MT FA = 142^{o}';
legend_cell{2}='MT FA = 426^{o}';
my_legend(legend_cell)

figure()

semilogx(protocol(1:2:end,2),abs(mt_measurement_F10perc_low(1:2:end)-mt_measurement_true(1:2:end))./mt_measurement_true(1:2:end).*100,'s','LineWidth',3)
hold on
semilogx(protocol(2:2:end,2),abs(mt_measurement_F10perc_low(2:2:end)-mt_measurement_true(2:2:end))./mt_measurement_true(2:2:end).*100,'rx','LineWidth',3)

my_xlabel('Offset Frequency (Hz)')
my_ylabel('% Difference vs. Reference')
axis([100 100000 -1 4])
title('Increase in F')
legend_cell{1}='MT FA = 142^{o}';
legend_cell{2}='MT FA = 426^{o}';
my_legend(legend_cell)

figure()
semilogx(protocol(1:2:end,2),abs(mt_measurement_kf10perc_low(1:2:end)-mt_measurement_true(1:2:end))./mt_measurement_true(1:2:end).*100,'s','LineWidth',3)
hold on
semilogx(protocol(2:2:end,2),abs(mt_measurement_kf10perc_low(2:2:end)-mt_measurement_true(2:2:end))./mt_measurement_true(2:2:end).*100,'rx','LineWidth',3)

my_xlabel('Offset Frequency (Hz)')
my_ylabel('% Difference vs. Reference')
axis([100 100000 -1 4])
title('Increase in kf')
legend_cell{1}='MT FA = 142^{o}';
legend_cell{2}='MT FA = 426^{o}';
my_legend(legend_cell)



%% End Program
%

%%% End timing %%%
howLong = toc;
time_result = ['Elapsed time is ', num2str((howLong-rem(howLong,60))/60), ' minutes and ' num2str(rem(howLong,60)) ' seconds'];
disp(time_result)
%%%