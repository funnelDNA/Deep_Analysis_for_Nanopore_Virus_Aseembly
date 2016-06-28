%% parameters to change
where2cut = 2.5; % number of sigmas
%%
hold off;
upLim = 0.01;
lowLim = 0;
binN = 150;
cutThickness = 1.5; % tentative cutoff, to make data cleaner before fitting. unit: sigma
%%
    tempDel = abs(aveDel);
    tempDel(tempDel<=lowLim | tempDel>=upLim) = [];
    tempDel = [tempDel; upLim; lowLim];
    
    cutP = 0.0037;
    [N P] = hist(tempDel,binN);
    
    N = N/sum(N);        
    Nt = N; Pt = P;

T3Position = 0;
T4Position = 0;

newT3P = 0.1;
newT4P = 0.1;
epslon = 0.0001;

mu = 0;     sigma = 0;
while abs(newT3P - T3Position) > epslon || abs(newT4P - T4Position) > epslon
    hold off;
    
    T3Position = newT3P;
    T4Position = newT4P;
    
    cutPL = mu - cutThickness*sigma;     cutPR = mu + cutThickness*sigma;
    N = Nt; P = Pt;
    N(P>cutPL & P<cutPR) = 0;   %P(P>cutPL & P<cutPR) = [];
    stdTempDel = std(tempDel)*2^0.5;
    coEff = fitHBV( N, P, stdTempDel );
    mu = coEff(2);
    sigma = coEff(3)/2^0.5;
    plot([mu - where2cut*sigma mu - where2cut*sigma],[0 coEff(1)],'k');
    plot([mu + where2cut*sigma mu + where2cut*sigma],[0 coEff(1)],'k');
    Pos1 = mu;    sigma1 = sigma;

    cutPL = mu - cutThickness*sigma;     cutPR = mu + cutThickness*sigma;
    N = Nt; P = Pt;
    N(P>cutPL & P<cutPR) = 0;   %P(P>cutPL & P<cutPR) = [];
    coEff = fitHBV( N, P, stdTempDel );
    mu = coEff(2);
    sigma = coEff(3)/2^0.5;
    plot([mu - where2cut*sigma mu - where2cut*sigma],[0 coEff(1)],'k');
    plot([mu + where2cut*sigma mu + where2cut*sigma],[0 coEff(1)],'k');
    Pos2 = mu;    sigma2 = sigma;

    newT3P = min([Pos1 Pos2]);
    newT4P = max([Pos1 Pos2]);
    
end
if Pos1<Pos2
    T3Position = Pos1;
    T3sigma = sigma1;
    T4Position = Pos2;
    T4sigma = sigma2;
else
    T3Position = Pos2;
    T3sigma = sigma2;
    T4Position = Pos1;
    T4sigma = sigma1;    
end

T3L = T3Position - where2cut*T3sigma;
T3R = T3Position + where2cut*T3sigma;
T4L = T4Position - where2cut*T4sigma;
T4R = T4Position + where2cut*T4sigma;

fprintf('T3 position: %f\t\t T4 position: %f\n',newT3P,newT4P);
fprintf('cutoff for sub T3: 0 - %f\n', T3Position - where2cut*T3sigma);
fprintf('cutoff for T3: %f - %f\n', T3Position - where2cut*T3sigma, T3Position + where2cut*T3sigma);
fprintf('cutoff for intermediates: %f - %f\n', T3Position + where2cut*T3sigma, T4Position - where2cut*T4sigma);
fprintf('cutoff for T4: %f - %f\n', T4Position - where2cut*T4sigma, T4Position + where2cut*T4sigma);

fprintf('Writing Excels...\n');

titles = {'aveAmp','DelI/I','aveDwell','pore2pore time','DelI/I_1st/2nd'};
nWrite = [fName '\' fName '_cutoffs.xls'];
xlswrite(nWrite,{'T3 Position'},1,'a1');        xlswrite(nWrite,T3Position,1,'b1');
xlswrite(nWrite,{'T4 Position'},1,'a2');        xlswrite(nWrite,T4Position,1,'b2');
xlswrite(nWrite,{'sub T3'},1,'a3');             xlswrite(nWrite,[0, T3L sum(Nt(Pt<T3L & Pt>=0))],1,'b3');
xlswrite(nWrite,{'T3'},1,'a4');                 xlswrite(nWrite,[T3L T3R sum(Nt(Pt<T3R & Pt>=T3L))],1,'b4');
xlswrite(nWrite,{'Intermediates'},1,'a5');      xlswrite(nWrite,[T3R T4L sum(Nt(Pt<T4L & Pt>=T3R))],1,'b5');
xlswrite(nWrite,{'T4'},1,'a6');                 xlswrite(nWrite,[T4L T4R sum(Nt(Pt<T4R & Pt>=T4L))],1,'b6');

fprintf('Done with %s.\n',fName);
fprintf('All Done.\n');