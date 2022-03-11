function outputTable = readTableIn(UniqCaptLR,UniqAllLR,TotCaptLR,TotAllLR,BadReads)
warning('off')
outputTable=table;
SmplNum = length(fieldnames(UniqAllLR.Chr1));

%% summing captured reads
outputTable{fieldnames(UniqCaptLR.GAL1),{'SUM_uniq'}} = zeros(SmplNum,1);
outputTable{fieldnames(TotCaptLR.GAL1),{'SUM_tot'}} = zeros(SmplNum,1);
outputTable{fieldnames(TotCaptLR.GAL1),{'SUM_LowMapQ'}} = zeros(SmplNum,1);
outputTable{fieldnames(TotCaptLR.GAL1),{'SUM_Unpair'}} = zeros(SmplNum,1);
outputTable{fieldnames(TotCaptLR.GAL1),{'SUM_FarPairs'}} = zeros(SmplNum,1);
for iGeneName=fieldnames(UniqCaptLR)'
    for iSample=fieldnames(UniqCaptLR.(iGeneName{1}))'
        % summing uniq reads
        outputTable{iSample{1},[iGeneName{1} '_uniq']} = 2 * length(UniqCaptLR.(iGeneName{1}).(iSample{1})); % *2 for paired reads
        outputTable{iSample{1},{'SUM_uniq'}} = outputTable{iSample{1},{'SUM_uniq'}} + outputTable{iSample{1},[iGeneName{1} '_uniq']};
        % summing total reads (including duplicates)
        outputTable{iSample{1},[iGeneName{1} '_tot']} = 2 * length(TotCaptLR.(iGeneName{1}).(iSample{1}));
        outputTable{iSample{1},{'SUM_tot'}} = outputTable{iSample{1},{'SUM_tot'}} + outputTable{iSample{1},[iGeneName{1} '_tot']};
        % summing badly mapped reads, unpaired, and far pairs
        outputTable{iSample{1},{'SUM_LowMapQ'}} = outputTable{iSample{1},{'SUM_LowMapQ'}} + length(BadReads.Capt.LowMapQ.(iGeneName{1}).(iSample{1}));
        outputTable{iSample{1},{'SUM_Unpair'}} = outputTable{iSample{1},{'SUM_Unpair'}} + length(BadReads.Capt.Unpair.(iGeneName{1}).(iSample{1}));
        outputTable{iSample{1},{'SUM_FarPairs'}} = outputTable{iSample{1},{'SUM_FarPairs'}} + length(BadReads.Capt.FarPairs.(iGeneName{1}).(iSample{1}));
    end
end
%% summing all reads (from all chromosomes)
outputTable{fieldnames(UniqAllLR.Chr1),{'SUMCHR_uniq'}} = zeros(SmplNum,1);
outputTable{fieldnames(TotAllLR.Chr1),{'SUMCHR_tot'}} = zeros(SmplNum,1);
outputTable{fieldnames(TotAllLR.Chr1),{'SUMCHR_LowMapQ'}} = zeros(SmplNum,1);
outputTable{fieldnames(TotAllLR.Chr1),{'SUMCHR_Unpair'}} = zeros(SmplNum,1);
outputTable{fieldnames(TotAllLR.Chr1),{'SUMCHR_FarPairs'}} = zeros(SmplNum,1);
for iChrName=fieldnames(UniqAllLR)'
    for iSample=fieldnames(UniqAllLR.(iChrName{1}))'
        % summing uniq reads
        outputTable{iSample{1},[iChrName{1} '_uniq']} = 2 * length(UniqAllLR.(iChrName{1}).(iSample{1})); % *2 for paired reads
        outputTable{iSample{1},{'SUMCHR_uniq'}} = outputTable{iSample{1},{'SUMCHR_uniq'}} + outputTable{iSample{1},[iChrName{1} '_uniq']};
        % summing total reads (including duplicates)
        outputTable{iSample{1},[iChrName{1} '_tot']} = 2 * length(TotAllLR.(iChrName{1}).(iSample{1}));
        outputTable{iSample{1},{'SUMCHR_tot'}} = outputTable{iSample{1},{'SUMCHR_tot'}} + outputTable{iSample{1},[iChrName{1} '_tot']};
        % summing badly mapped reads, unpaired, and far pairs
        outputTable{iSample{1},{'SUMCHR_LowMapQ'}} = outputTable{iSample{1},{'SUMCHR_LowMapQ'}} + length(BadReads.All.LowMapQ.(iChrName{1}).(iSample{1}));
        outputTable{iSample{1},{'SUMCHR_Unpair'}} = outputTable{iSample{1},{'SUMCHR_Unpair'}} + length(BadReads.All.Unpair.(iChrName{1}).(iSample{1}));
        outputTable{iSample{1},{'SUMCHR_FarPairs'}} = outputTable{iSample{1},{'SUMCHR_FarPairs'}} + length(BadReads.All.FarPairs.(iChrName{1}).(iSample{1}));
    end
end
warning('on')