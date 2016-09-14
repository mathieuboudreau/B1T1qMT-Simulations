function mxya = mtspgr_rp2_model(p, sat_angles, sat_offsets, sat_pulses, lineshape, excite_angle, tr)
%
% function mxya = mtspgr_rp2_model(p, sat_angles, sat_offsets, sat_pulses, lineshape, excite_angle, tr)
%
%	calculate resulting signal from a pulsed MTSPGR sat experiment according to Sled & Pike's RP model
%       (and compare to John's implementation if desired - hardcoded)
%
% inputs:
%	sat_angles:	sat pulse flip angles (degrees)
%	sat_offsets:    sat pulse frequency offsets (Hz)
%	excite_angles:	excitation pulse flip angle (degrees)
%	p: 	    input parameters (1 = liquid/lorentzian, 2 = solid/specify)
%		    [f r t2a t2b r1a r1b]
%
%

p

ra = p(5);
t2a = p(3);
m0a = 1;

rb = p(6);
t2b = p(4);
m0b = p(1);

m0 = [m0a m0b];

kba = p(2);
kab = kba*m0b/m0a;

%lineshape_fcn = [lineshape '_lineshape'];
lineshape_fcn = lineshape;

mza = zeros(length(sat_angles), length(sat_offsets));
mxya = zeros(length(sat_angles), length(sat_offsets));

s_mt = zeros(length(sat_angles),length(sat_offsets));
s_excite = cos(excite_angle*pi/180);



%---- compute direct effect

for del = 1:length(sat_offsets)
  for fl = 1:length(sat_angles)
  end
end



for fl = 1:length(sat_angles)
  for del = 1:length(sat_offsets)
    sf = sim_sf_o(sat_pulses{fl,del},1/ra,t2a,sat_offsets(del));
    s_mt(fl,del) = sf;

    % John's method computes the saturation rate based on CW approximation, for power spread out across entire TR - not pulse duration!
    wr(fl,del) = cw_rfrate(sat_pulses{fl,del}, t2b, sat_offsets(del), lineshape_fcn);	

    % this is to reverse the scaling that happens within the computation for wr (see above)!  bad implementation?
    tau = pulse_duration(sat_pulses{fl,del});

    we = wr(fl,del)*tr/tau;
    
    % the matrix way
    s1= diag([sf 1]);
    s2 = diag([s_excite 1]);
        
    mss = get_steady_state(m0,ra,rb,kab,kba,we);

    afp = [-(ra + kab), kba; kab, -(rb + kba)];
    expafpt = expm(afp*(tr-tau));

    acw = [-(ra + kab), kba; kab, -(rb + kba + we)];
    expacwt = expm(acw*tau/2);

        
    %expacwt
    %expafpt
    %s1
    %s2
    %mss'
    %m0'
    
    m = (inv(eye(2) - expacwt*expafpt*expacwt*s2*s1))*((expacwt*expafpt - expacwt*expafpt*expacwt + eye(2) -expacwt)*mss' + expacwt*(eye(2)-expafpt)*m0');
    mza(fl,del) = m(1);
    mxya(fl,del) = mza(fl,del)*sin(excite_angle*pi/180);
        
    %--- using Sled's code
    %[a b] = mtspgr_rp2(excite_angle, kab, m0b, tr, [1/ra 1/rb], [s_mt(fl,del) 1] , wr(fl,del), t2b, sat_offsets(del)*2*pi, tau);
    %mz_sled(fl,del) = a(1,:);
    %mxy_sled(fl,del) = b;
  end 
end


return

%--- "no sat" baseline data

% MT pulse amplitude = 0
omega_1 = 0;
we = 0;

% the matrix way
s1= diag([1 1]);
s2 = diag([s_excite 1]);
        
mss = get_steady_state(m0,ra,rb,kab,kba,we);

afp = [-(ra + kab), kba; kab, -(rb + kba)];
expafpt = expm(afp*(tr-tau));

acw = [-(ra + kab), kba; kab, -(rb + kba + we)];
expacwt = expm(acw*tau/2);
         
m = (inv(eye(2) - expacwt*expafpt*expacwt*s2*s1))*((expacwt*expafpt - expacwt*expafpt*expacwt + eye(2) -expacwt)*mss' + expacwt*(eye(2)-expafpt)*m0');
mza_nosat = m(1);
mxya_nosat = mza_nosat*sin(excite_angle*pi/180)*sf;

%--- using Sled's code
[a b] = mtspgr_rp2(excite_angle(expt), kab, m0b, tr, [1/ra 1/rb], [1 1] , 0, t2b, sat_offsets(del)*2*pi, tau);
mz_sled_nosat = a(1,:);
mxy_sled_nosat = b;




%-----------------------------------------------------------------------------------------------------
function mss = get_steady_state(m0,r1f,r1r,kf,kr,w)

denom = r1r*r1f + r1r*kf + r1f*kr + w*r1f + w*kf;

mss(1) = m0(1)*(r1r*kf + r1r*r1f + r1f*kr + w*r1f)/denom;
mss(2) = m0(2)*(r1r*r1f+ r1r*kf + r1f*kr)/denom;
