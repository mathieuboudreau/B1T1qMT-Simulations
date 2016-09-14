function [dM, out2, out3] = bloch_1pool(t, M, flag, sequence, T1, T2, dw)
%
%  function dM = bloch_1pool(t, M, flag, sequence, T1, T2, Ka, F, dw)
%
%  ODE file implementing equations of motion for single pool
%  governed by the Bloch equations.  B1 waveform is transformed to
%  take into account rotating frame.
%
%  t          : time in seconds
%  M          : column vector with three magnetization components
%               [MxA MyA MzA]'
%  sequence   : sequence object such that omega1(sequence,t) = gamma*B1(t)
%  T1         : longitudinal relaxation
%  T2         : transverse relaxation
%  dw         : offset frequency of rotating frame = w - w0 (radians)

if nargin < 3 | isempty(flag)

  if(isa(sequence, 'double'))
    w1 = 0;
  else
    w1 = exp(sqrt(-1)*t*dw)*omega1(sequence,t);
  end
  w1x = real(w1);
  w1y = imag(w1);
  dM = zeros(3,size(M,2));

  MxA = M(1,:);
  MyA = M(2,:);
  MzA = M(3,:);

  dM(1,:) = -MxA/T2(1) - dw*MyA - w1y*MzA;
  dM(2,:) = -MyA/T2(1) + dw*MxA + w1x*MzA;
  dM(3,:) = (1 - MzA)/T1(1) + w1y*MxA - w1x*MyA;
%  dM(3,:) = w1y*MxA - w1x*MyA;
  
else
  switch(flag)
    case 'init'                           % Return default [tspan,y0,options].
      dM = [0 1];
      out2 = [0 0 1];
      out3 = [];
    case 'jacobian'
      if(isa(sequence, 'double'))
        w1 = 0;
      else
        w1 = exp(sqrt(-1)*t*dw)*omega1(sequence,t);
      end
      w1x = real(w1);
      w1y = imag(w1);
      dM = zeros(3,3);
      dM(1,:) = [ -1/T2(1) -dw       -w1y ];
      dM(2,:) = [    dw    -1/T2(1)   w1x ];
      dM(3,:) = [    w1y    -w1x  -1/T1(1)];
   
    otherwise
      error(['Unknown flag ''' flag '''.']);
  end
end

