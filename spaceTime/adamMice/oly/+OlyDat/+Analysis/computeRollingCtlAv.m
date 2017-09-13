function [data tvec ifo] = computeRollingCtlAv(data,statFld,timeFld,nTime)
% data = computeRollingCtlAv(data,statFlds,timeFld,nTime)
% Compute rolling average of stats over time.
% timeFld: time-related field, eg 'day_num', 'week_num', etc. data.timeFld
% should be an integer.
% nTime: number of time units to average over.

fprintf(1,'Computing rolling averages...\n');

tfCtl = OlyDat.Analysis.getNumericScalarFieldSafe(data,'tfCtl');

allTime = OlyDat.Analysis.getNumericScalarFieldSafe(data,timeFld);
maxTime = max(allTime);
tvec = 1:maxTime;

for t = tvec
    if mod(nTime,2)==0
        tWin = t-nTime/2:t+nTime/2-1;
    else
        tWin = t-(nTime-1)/2:t+(nTime-1)/2;
    end
    
    tfThisWin = ismember(allTime,tWin);
    yCtlThisWin = OlyDat.Analysis.getNumericScalarFieldSafe(data(tfCtl&tfThisWin),statFld);
    yGmrThisWin = OlyDat.Analysis.getNumericScalarFieldSafe(data(~tfCtl&tfThisWin),statFld);
    if isempty(yCtlThisWin)
        yCtlThisWin = yCtlThisWin(:);
    end
    if isempty(yGmrThisWin)
        yGmrThisWin = yGmrThisWin(:);
    end        
    assert(isequal(size(yCtlThisWin),[nnz(tfCtl&tfThisWin) 1]));
    assert(isequal(size(yGmrThisWin),[nnz(~tfCtl&tfThisWin) 1]));
    
    ifo(t,1).win = tWin; %#ok<*AGROW>
    ifo(t,1).NCtl = numel(yCtlThisWin);
    ifo(t,1).meanCtl = nanmean(yCtlThisWin);
    ifo(t,1).mdnCtl = nanmedian(yCtlThisWin);
    ifo(t,1).stdCtl = nanstd(yCtlThisWin);
    [~,~,ifo(t,1).meanCtlCI(1,:),ifo(t,1).sdCtlCI(1,:)] = normfit(yCtlThisWin(~isnan(yCtlThisWin)));
    ifo(t,1).NGmr = numel(yGmrThisWin);
    ifo(t,1).meanGmr = nanmean(yGmrThisWin);
    ifo(t,1).mdnGmr = nanmedian(yGmrThisWin);
    ifo(t,1).stdGmr = nanstd(yGmrThisWin);
    
    avFldBase = sprintf('%s__%d%s__',statFld,nTime,timeFld);
    tfThisT = allTime==t;
    if nnz(tfThisT)>0
        [data(tfThisT).([avFldBase 'tWin'])] = deal(mat2str(ifo(t).win));
        [data(tfThisT).([avFldBase 'NCtl'])] = deal(ifo(t).NCtl);
        [data(tfThisT).([avFldBase 'meanCtl'])] = deal(ifo(t).meanCtl);
        [data(tfThisT).([avFldBase 'mdnCtl'])] = deal(ifo(t).mdnCtl);
        [data(tfThisT).([avFldBase 'stdCtl'])] = deal(ifo(t).stdCtl);
        
        [data(tfThisT).([avFldBase 'NGmr'])] = deal(ifo(t).NGmr);
        [data(tfThisT).([avFldBase 'meanGmr'])] = deal(ifo(t).meanGmr);
        [data(tfThisT).([avFldBase 'mdnGmr'])] = deal(ifo(t).mdnGmr);
        [data(tfThisT).([avFldBase 'stdGmr'])] = deal(ifo(t).stdGmr);        
    end
end
end

