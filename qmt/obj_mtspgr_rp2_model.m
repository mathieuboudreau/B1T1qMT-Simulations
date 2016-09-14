function z = obj_mtspgr_rp2_model(params, measurements, t1obs, lineshape, sat_angles, offsets, sat_pulses, excite_angle, tr, r1b, params_expnt, dim_sat_angles, protocolFlag )
%
% z = obj_mtspgr_rp2_model(params, measurements, t1obs, lineshape, sat_angles, offsets, sat_pulses, excite_angle, r1b)
%
% inputs:
%  params(1) = f
%  params(2) = r
%  params(3) = t2a
%  params(4) = t2b
%  params(7) = scale
%
%  r1a computed from t1obs and t1b
%
%   objective function for fitting the RP model of pulsed-sat-SPGR
%
%
%

%Rescale parameters
params = params.*10.^(params_expnt); 

% compute r1a, assuming r1b = 1
params(6) = r1b;
rd = params(6) - 1/t1obs;
kf = params(2)*params(1);
params(5) = 1/t1obs  - kf * rd /(rd + params(2));


%Memory pre-allocation
estimate_mt = zeros(dim_sat_angles, length(offsets), length(tr));
estimate = zeros(dim_sat_angles, length(offsets), length(tr));

% Pre-allocate memory.
switch protocolFlag
    case 'sled' 
        estimate_mt = zeros(dim_sat_angles, length(offsets), length(tr));
        estimate = zeros(dim_sat_angles, length(offsets), length(tr));
    case 'custom'
        estimate_mt = zeros(dim_sat_angles, length(offsets), length(tr));
        estimate = cell(length(tr),1);
end


for expt = 1:length(tr)
    switch protocolFlag
        case 'sled'
            estimate_mt(:,:,expt) = mtspgr_rp2_model(params, sat_angles(expt,:), offsets, sat_pulses(:,:,expt), lineshape, excite_angle(expt), tr(expt)); 
            estimate_baseline = max(max(estimate_mt(:,:,expt))); 
            estimate(:,:,expt) = estimate_mt(:,:,expt)/estimate_baseline; 
        case 'custom' 
            custom_sat_angles = cell2mat(sat_angles{expt});
            estimate_baseline = mtspgr_rp2_model_baseline(params, custom_sat_angles, offsets, sat_pulses(:,:,expt), lineshape, excite_angle(expt), tr(expt));
            estimate_mt = mtspgr_rp2_model(params, custom_sat_angles, offsets, sat_pulses(:,:,expt), lineshape, excite_angle(expt), tr(expt));

            estimate{expt} = {estimate_mt/estimate_baseline};
    end
end

switch protocolFlag
    case 'sled'
        z = sum(sum(sum((measurements - estimate).^2)));
    case 'custom'
        for  expt = 1:length(tr)
            z(expt) = sum(sum((cell2mat(measurements{expt}) - cell2mat(estimate{expt})).^2));
            
        end
        z=sum(z);
end
        
