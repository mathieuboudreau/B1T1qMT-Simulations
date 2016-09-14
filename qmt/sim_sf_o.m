function Sf = sim_sf_o(pulse, T1, T2, offset, M0)
%
%  function Sf = sim_sf(pulse, T1, T2, offset)
%
%  compute direct effect on free water pool 
%
%  offset      :  offset frequency of rotating frame used for simulation
%  M0          :  current Z magnetization

if(nargin == 4)
  M0 = [0 0 1.0]';
else
  M0 = [0 0 M0]';
end

dw = -2*pi*offset;
tau = pulse_duration(pulse);

%options = odeset('Vectorized', 'on', 'Jacobian', 'on', 'AbsTol', 1e-6, 'RelTol', 1e-9, 'MaxOrder', 2);
options = odeset('Vectorized', 'on', 'Jacobian', 'on', 'AbsTol', 1e-6, 'RelTol', 1e-9, 'MaxOrder', 2);
%

if(abs(offset) < 5e3)
  [T,M] = ode45('bloch_1pool', [0 tau], M0, options, pulse, T1, T2, dw);
elseif(abs(offset) < 1e4)
  [T,M] = ode45('bloch_1pool', [0 tau], M0, options, pulse, T1, T2, dw);
  [T_,M_] = ode15s('bloch_1pool', [0 tau], M0, options, pulse, T1, T2, dw);
  a = (log10(offset)-log10(5e3))/(log10(5e4)-log10(5e3));
  M = (1-a)*M(length(T),:) + a*M_(length(T_),:);
else
  [T,M] = ode15s('bloch_1pool', [0 tau], M0, options, pulse, T1, T2, dw);
end

Sf = M(size(M,1),3);
return
figure(1);
hold off;
plot(T, M(:,3));

figure(2);
hold off;
plot(T, M(:,1));

figure(3);
hold off;
plot(T, M(:,2));
