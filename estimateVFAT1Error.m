function [T1error, measuredSignal, Sy, Sx] = estimateVFAT1Error(T1,TR,FAs, faErrorRange)
%ESTIMATEVFAT1ERROR Summary of this function goes here
%   T1: scalar in ms
%   TR: scalar in ms
%   FAs: array in deg
%   faErrorRange: array in relative amplitude

M0 = 10;
measuredSignal = vfaSignal(T1,TR,FAs,M0);

%kFactor = 1-exp(-TR/T1);
%T1error = TR./( log(cosd(FA.*faErrorRange)) - log((1 - (sind(FA.*faErrorRange)./measuredSignal ) .* kFactor .* M0) ));

for ii=1:length(faErrorRange)
    Sy = measuredSignal./sind(FAs.*faErrorRange(ii));
    Sx = measuredSignal./tand(FAs.*faErrorRange(ii));

    p=polyfit(Sx,Sy,1);

    T1error(ii) = -TR./log(p(1));
end

end

function [signal] = vfaSignal(T1,TR,FA,M0)
% T1 in ms
% TR in ms
% FA in deg

    E1 = exp(-TR/T1);
    signal=M0 .* (1-E1)./(1-cosd(FA).*E1) .* sind(FA);

end