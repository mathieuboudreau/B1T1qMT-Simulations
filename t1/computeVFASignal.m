function [ VFASignal ] = computeVFASignal( seqParam )
%COMPUTEVFASIGNAL Summary of this function goes here
%   Detailed explanation goes here
    
    TR = seqParam.TR;
    T1 = seqParam.T1;
    FlipAngles = seqParam.FlipAngles;
    equilibMagn =1;
    
    VFASignal = equilibMagn.*((1-exp(-TR./T1)).*sind(FlipAngles))./(1-exp(-TR./T1).*cosd(FlipAngles)); 
    
end

