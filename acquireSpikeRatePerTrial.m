function sr = acquireSpikeRatePerTrial(st,timeBin,tAxis)
    % convert spike timing data into spike rate for each trial
    
    sr = nan(numel(st),numel(tAxis)); 
    
    for ntrial = 1:numel(st)
        spkTime      = st{ntrial};
        sr(ntrial,:) = hist(spkTime,tAxis)/timeBin; 
    end
    
end

