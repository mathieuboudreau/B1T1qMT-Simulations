function p = fit_mtspgr_rp2_model(sat_angles, sat_offsets, sat_duration, excite_angle, lineshape, measurements, t1obs, tr, r1_solid, protocolFlag)
%
% p = fit_mtspgr_rp2_model(sat_angles, sat_offsets, sat_duration, excite_angle, lineshape, measurements, t1obs, tr, r1_solid)
%
%	A function to fit pulsed saturation data with the Sled RP-approximation analytical model
%		and return two-pool model parameters
%
% output
%	p = structure with fit parameters
%
% input
%	saturation angles and offsets
%	solid pool lineshape to be used
%	signal measurements
%	value of t1obs measure independently
%	[fixed value of R1r to use]
%
%------------------------------------------------------------------------------


% create initial guess
p0.f = 0.122;  
p0.t2 = [.0272 12e-6];
p0.r = 32.541;


if ~exist('r1_solid','var')
    p0.r1(2) = 1;
else
    p0.r1(2) = r1_solid;
end

% construct initial guess for r1a, assuming r1b = 1
%p0.r1(2) = r1_solid;
rd = p0.r1(2) - 1/t1obs;
kf = p0.r*p0.f;
p0.r1(1) = 1/t1obs  - kf * rd /(rd + p0.r);

guess = [p0.f p0.r p0.t2(1) p0.t2(2)];
[guess_mant guess_expnt] = mantexpnt(guess);

% initialize fitting options structure
myopts = optimset('fminsearch');
%optimget(myopts, 'MaxFunEvals')
myopts.MaxFunEvals = 200*5;

myopts.TolX = 0.1; 
myopts.Display = 'iter'; 

lb = [0   0   0   0   0   0   0];
ub = [1 Inf Inf Inf Inf Inf Inf];

% Pre-allocate memory.
my_pulses = cell(length(sat_angles),length(sat_offsets),length(tr));

% Pulse setup
switch protocolFlag
    case 'custom'
        for forIndex = 1:length(sat_angles)

            temp_sat_angles = sat_angles{forIndex};
            temp_sat_angles = cell2mat(temp_sat_angles);
            each_dim_sat_angles(forIndex) = length(temp_sat_angles);
            clear temp_sat_angles
        end
    dim_sat_angles= sum(each_dim_sat_angles);
    
end

for expt = 1:length(tr)
  for del = 1:length(sat_offsets)
    
    switch protocolFlag
        case 'sled'
            fl_max = length(sat_angles);
        case 'custom'
            fl_max = each_dim_sat_angles(expt);
    end
            
    for fl = 1:fl_max
       switch protocolFlag
           case 'sled' 
                my_pulses{fl,del,expt} = gaussian_hann(sat_angles(expt,fl),sat_duration(expt), sat_offsets(del),tr(expt),200); 
           case 'custom'
                temp_sat_angles = cell2mat(sat_angles{expt});
                my_pulses{fl,del,expt} = gaussian_hann(temp_sat_angles(fl),sat_duration(expt), sat_offsets(del),tr(expt),200); 
                clear temp_sat_angles;
       end
    end
  end
end



'pulses done, now fit'

%[result, resnorm, res, flag] = lsqnonlin(@obj_cw_model, guess, lb, ub,[],
%measurements, t1obs, lineshape, amplitudes, sat_offsets);
[result, fval, flag] = fminsearch(@obj_mtspgr_rp2_model, guess_mant, myopts, measurements, t1obs, lineshape, sat_angles, sat_offsets, my_pulses, excite_angle, tr, p0.r1(2),guess_expnt, dim_sat_angles, protocolFlag); 
result = result.*(10.^guess_expnt); 


%result
%fval
%flag

%while ~flag
%    % increase the tolerance
%    myopts = optimset(myopts,'MaxFunEvals', myopts.MaxFunEvals*2)
%    [result, resnorm, res, flag] = lsqnonlin(@obj_cw_model, guess, lb, ub, myopts, measurements, t1obs, lineshape, amplitudes, sat_offsets);
%    flag
%end



while ~flag
    % increase the tolerance
    myopts = optimset(myopts,'MaxFunEvals', myopts.MaxFunEvals*2);
    flag
    [result, fval, flag] = fminsearch(@obj_mtspgr_rp2_model, guess, myopts, measurements, t1obs, lineshape, sat_angles, sat_offsets, my_pulses, excite_angle, tr, p0.r1(2));
end


p = p0;

p.f = result(1);
p.r = result(2);

rd = p.r1(2) - 1/t1obs;
kf = p.r*p.f;
p.r1(1) = 1/t1obs  - kf * rd /(rd + p.r);

p.t2(1) = result(3);
p.t2(2) = result(4);
%p.scale = result(7);


