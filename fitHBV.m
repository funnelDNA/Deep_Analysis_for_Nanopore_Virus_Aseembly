function [ coEff ] = fitHBV( N, P, stdTempDel )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    [maxN maxNP] = max(N);
    [f gof] = fit(P',N','gauss1','startpoint',[maxN P(maxNP) stdTempDel],'lower',[0 0 0.0],'upper',[1 0.02 0.01]);
    coEff = coeffvalues(f);    
    mu = coEff(2);
    sigma = coEff(3)/2^0.5;
    plot(f,'r'); hold on;
    plot(P,N,'o');
%    fprintf('T = 4 position: %f\t\t',mu);
    c = sigma/2^0.5;
    area = coEff(1)*sigma*3.1415926^0.5/(1/150*0.01)*2^0.5;
end

