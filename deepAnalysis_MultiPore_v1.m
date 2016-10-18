clear;
allFiles = dir('*.abf');
cutE = 3;
nPores = 3;
for ii = 1:length(allFiles)
    fName = allFiles(ii).name;
    fprintf('Analyzing %s...\n',fName);
    fName(end-3:end) = [];
    load ([fName '\ExperimentData.mat']);
    startstop = cell2mat(transpose([{EventDatabase.StartAndEndPoint}]));
    amp = cell2mat(transpose([{ EventDatabase.AllLevelFits}]));
    baseline = cell2mat(transpose([{ EventDatabase.deltai}]));
    start = startstop(:,1)/SamplingFrequency;
    stop  = startstop(:,2)/SamplingFrequency;
    ppTime = start(2:end) - stop(1:end-1);
    ppTimeBak = ppTime;
    upLim = 0.1;
    lowLim = 0;
    for jj = 1:2
        ppTime = ppTimeBak;
        ppTime(ppTime>upLim) = [];
        ppTime(ppTime<=lowLim) = [];
        [N P] = hist(log(ppTime),100);
        
        f = fit(P',N','gauss1','startpoint',[1 log(median(ppTime)) log(mean(ppTime))]);
        fitY = f.a1*exp(-((P'-f.b1)/f.c1).^2);
        coEff = coeffvalues(f);
        mu = coEff(2);
        sigma = coEff(3)/2^0.5;
        rightP = min((mu+cutE*sigma));
        if jj == 1
            upLim = exp(rightP)*2;
        else
            upLim = 0;
        end
    end
    %%
    plot(exp(P'),fitY,exp(P'),N','o');
    hold on;
    %bar(exp(P'),N');
    
    leftP  = max((mu-cutE*sigma));
    rightP = min((mu+cutE*sigma));
    plot([exp(leftP) exp(leftP)],[0 max(N)/2],'color','r');
    plot([exp(rightP) exp(rightP)],[0 max(N)/2],'color','r');
    xlabel('pore to pore time(s)','fontsize',20);
    ylabel('counts','fontsize',20);
    %set(gca,'ylim', [-1 coEff(1)*2],'xlim',[0 mu*4],'fontsize',20);
    hold off;
    %%
    leftCut = exp(leftP);
    rightCut = exp(rightP);
    fprintf('Left cutoff = %f , right cutoff = %f\n',leftCut, rightCut);
    isManualCutoff = input('Do you want to set cutoffs manually?(Y/N)\n','s');
    while (isManualCutoff ~= 'Y') && isManualCutoff ~= 'y' && isManualCutoff ~= 'N' && isManualCutoff ~= 'n'
        isManualCutoff = input('Please type Y or N\n','s');
    end
    if isManualCutoff == 'Y' || isManualCutoff == 'y'
        leftCut = input('please type in Left cutoff...\n');
        rightCut = input('please type in right cutoff...\n');        
    end
    %%
    ppTime = ppTimeBak;
    dataIndex = find(ppTime<=rightCut & ppTime>=leftCut);
    aveAmp = zeros(length(dataIndex),1);
    aveDel = aveAmp;
    aveDwell = aveAmp;
    p2p    = aveAmp;
    P1overP2 = aveAmp;
    del = amp./baseline;
    dwell = stop - start;
    jj = 1;
    ll = 1;
    while (jj<length(dataIndex))
        kk = 1;
%        if (jj+nPores >= length(dataIndex))
%            break;
%        end
        
        tempAmp = zeros(1,nPores);
        tempDel = zeros(1,nPores);
        tempDwell = zeros(1,nPores);
        tempP2p = zeros(1,nPores-1);
        endEarly = 0;
        for kk = 1:nPores-1
            if jj > length(dataIndex) continue; end
            tempAmp(kk) = amp(dataIndex(jj));
            tempDel(kk) = del(dataIndex(jj));
            tempDwell(kk) = dwell(dataIndex(jj));
            tempP2p(kk) = ppTime(dataIndex(jj));
            if (kk<nPores-1) && jj < length(dataIndex)
                if dataIndex(jj+1) > dataIndex(jj)+1
                    endEarly = 1;
                    jj = jj+1;
                    break;
                end
            end
            jj = jj + 1;
        end
        if endEarly == 1
            continue;
        end
        if (jj<=length(dataIndex) && dataIndex(jj)-dataIndex(jj-1)<=1)
            jj = jj + 1;
            while (jj<=length(dataIndex) && dataIndex(jj) == dataIndex(jj-1)+1)
                jj = jj+1;
            end
            continue;
        end
        if endEarly == 0
            tempAmp(end) = amp(dataIndex(jj-1)+1);
            tempDel(end) = del(dataIndex(jj-1)+1);
            tempDwell(end) = dwell(dataIndex(jj-1)+1);
            %tempP2p(end) = ppTime(dataIndex(jj-1));
            
            aveAmp(ll) = sum(tempAmp)/nPores;
            aveDel(ll) = sum(tempDel)/nPores;
            aveDwell(ll) = sum(tempDwell)/nPores;
            p2p(ll) = sum(tempP2p)/(nPores-1);
            ll = ll + 1;
        end
        %P1overP2(jj) = del(dataIndex(jj))/del(dataIndex(jj)+1);
    end
    
    aveAmp(ll:end) = [];
    aveDel(ll:end) = [];
    aveDwell(ll:end) = [];
    p2p(ll:end) = [];
    P1overP2(ll:end) = [];
    
    fprintf('Writing Excels...\n');

    titles = {'aveAmp','DelI/I','aveDwell','pore2pore time','DelI/I_1st/2nd'};
    nWrite = [fName '\' fName '_new.xls'];
    xlswrite(nWrite,titles,1,'a1');
    xlswrite(nWrite,[abs(aveAmp) abs(aveDel) aveDwell p2p P1overP2],1,'a2');
    xlswrite(nWrite,{'','p-p time(s)';'High',rightCut;'Low',leftCut},1,'g2:h4');
    xlswrite(nWrite,{'counts',length(aveAmp);'counts/min',length(aveAmp)/(stop(end)/60);'total count',length(amp)/2},1,'g6:h8');
    %xlswrite(nWrite,{'ave P1/P2',mean(P1overP2);'std P1/P2',std(P1overP2)},1,'g10:h11');
    fprintf('Done with %s.\n',fName);
end
fprintf('All Done.\n');