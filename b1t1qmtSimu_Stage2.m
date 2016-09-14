function [] = b1t1qmtSimu_Stage2(outputDir,mt_measurement,vol_dims,vol_resolution,protocol, sled_trs, sled_excite_flips)
%B1T1QMTSIMU_STAGE2 Creates Minc files for the qMT fitting to be performed, according to the B1 and T1 maps defined within this file. 
%
%   -The qMT images were are flat (every pixel the same value) and equal
%   to the simulated data from main_qmt_simu.m.
%
%   -The B0 measurements are flat and equal to zero.
%
%   -The B1 map is a 1D gradient along the left-right direction.
%
%      B1_min           B1_max
%           _ _ _ _ _ _ _ _ 
%          |               |
%          |               |
%          |               |
%          |               |
%          |               |
%          |               |
%          |               |
%          |_ _ _ _ _ _ _ _|
%    
%
%   -The T1 map is fitted from simulated VFA data for a range of T1 1D
%   gradient along the up-down direction. The data is fitted using the B1
%   map above.
%
%
%      (B1_min)         (B1_max)
%           _ _ _ _ _ _ _ _ 
%  T1_min  |               |
%          |               |
%          |               |
%          |               |
%          |               |
%          |               |
%          |               |
%  T1_max  |_ _ _ _ _ _ _ _|
%
%
%   -The dR1 images are required by the qMT fitting software, but are set 
%   to 0.
%
%   ***NOTE THAT NIAK REORIENTS THE IMAGES, SO THAT THE MINC FILES WILL NOT
%   HAVE THE GRADIENTS IN THE SAME DIRECTION***
%
%   Author: Mathieu Boudreau
%   Date Created: June 2014
%   Date Last Modified: June 20th 2014.
%

%% Create mask
%

mask=ones(vol_dims);
MASK = minc;
MASK.fileName=[outputDir 'mask.mnc'];
MASK=MASK.setVolume(mask);
MASK=MASK.setResolution(vol_resolution);
MASK.saveMinc;

%% Create qMT Minc volumes
%
%Set No-MT volumes
qmt_m0_vol=ones(vol_dims);

if length(sled_trs)==1 %Use for UK protocol
    qMT_MINC_M0=minc;
    qMT_MINC_M0.fileName=[outputDir 'qMT/qmt_0.mnc'];
    qMT_MINC_M0=qMT_MINC_M0.setVolume(qmt_m0_vol);
    qMT_MINC_M0=qMT_MINC_M0.setResolution(vol_resolution);
    qMT_MINC_M0=qMT_MINC_M0.setTR(sled_trs);
    qMT_MINC_M0=qMT_MINC_M0.setFA(sled_excite_flips);
    qMT_MINC_M0.saveMinc;
else
    for ii=1:length(sled_trs)
    qMT_MINC_M0=minc;
    qMT_MINC_M0.fileName=[outputDir 'qMT/qmt_0_tr_' num2str(sled_trs(ii)) '.mnc'];
    qMT_MINC_M0=qMT_MINC_M0.setVolume(qmt_m0_vol);
    qMT_MINC_M0=qMT_MINC_M0.setResolution(vol_resolution);
    qMT_MINC_M0=qMT_MINC_M0.setTR(sled_trs(ii));
    qMT_MINC_M0=qMT_MINC_M0.setFA(sled_excite_flips(ii));
    qMT_MINC_M0.saveMinc;                
    end
end


% Set MT weighted volumes
for qmtCount=1:length(mt_measurement)
   qMT_MINC(qmtCount)=minc;
   qMT_MINC(qmtCount).fileName=[outputDir 'qMT/qmt_' num2str(qmtCount) '.mnc'];
   qMT_MINC(qmtCount)=qMT_MINC(1,qmtCount).setVolume(mt_measurement(qmtCount).*ones(vol_dims));
   qMT_MINC(qmtCount)=qMT_MINC(1,qmtCount).setResolution(vol_resolution);
   qMT_MINC(qmtCount)=qMT_MINC(1,qmtCount).setTR(protocol(qmtCount,4));
   qMT_MINC(qmtCount)=qMT_MINC(1,qmtCount).setFA(protocol(qmtCount,5));
   qMT_MINC(qmtCount).saveMinc;
end

%% Create B0 Minc volume
%

b0_vol=zeros(vol_dims);

B0 = minc;
B0.fileName=[outputDir 'B0/b0.mnc'];
B0=B0.setVolume(b0_vol);
B0=B0.setResolution(vol_resolution);
B0.saveMinc;

%% Create B1 Minc volume
%

b1_range=[0.5 1.5];
%b1_range=[1 1];

b1_vol=repmat(linspace(b1_range(1),b1_range(2),vol_dims(1)),vol_dims(2),1);

B1 = minc;
B1.fileName=[outputDir 'B1/b1.mnc'];
B1=B1.setVolume(b1_vol);
B1=B1.setResolution(vol_resolution);
B1.saveMinc;

%% Create T1 Minc volume
%

%CREATE VFA ACQ FILES
t1_range=linspace(0.5,1.3,vol_dims(1));

seqParam.TR=0.015;
seqParam.FlipAngles=[3 20];

for ii=1:vol_dims(1)
    seqParam.T1=t1_range(ii);
    [ VFASignal(ii,:) ] = computeVFASignal( seqParam );
end

%SAVE FIRST FA VFA IMAGE
vfa_vol{1}=repmat(VFASignal(:,1),1,vol_dims(1));
VFA1 = minc;
VFA1.fileName=[outputDir 'T1/vfa_1.mnc'];
VFA1=VFA1.setVolume(vfa_vol{1});
VFA1=VFA1.setResolution(vol_resolution);
VFA1=VFA1.setTR(seqParam.TR);
VFA1=VFA1.setFA(seqParam.FlipAngles(1));
VFA1.saveMinc;


%SAVE SECOND FA VFA IMAGE
vfa_vol{2}=repmat(VFASignal(:,2),1,vol_dims(1));
VFA2 = minc;
VFA2.fileName=[outputDir 'T1/vfa_2.mnc'];
VFA2=VFA2.setVolume(vfa_vol{2});
VFA2=VFA2.setResolution(vol_resolution);
VFA2=VFA2.setTR(seqParam.TR);
VFA2=VFA2.setFA(seqParam.FlipAngles(2));
VFA2.saveMinc;

% Fits T1 maps, saves Minc files.
qt1_vfa_lin_fit({[outputDir 'T1/vfa_1.mnc'],[outputDir 'T1/vfa_2.mnc']}, seqParam.FlipAngles, seqParam.TR, [outputDir 'T1/t1.mnc'], [outputDir 'T1/m0.mnc'], [outputDir 'B1/b1.mnc'], [outputDir 'mask.mnc'])


%% Create dR1 Minc volume
%

dr1_vol=zeros(vol_dims);

dR1 = minc;
dR1.fileName=[outputDir 'T1/dr1.mnc'];
dR1=dR1.setVolume(dr1_vol);
dR1=dR1.setResolution(vol_resolution);
dR1.saveMinc;

end

