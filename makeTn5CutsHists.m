function [HistUniqNormLoci,HistUniqNormChr,HistTotNormLoci,HistTotNormChr,StartStop,FragLenHist] ...
    = makeTn5CutsHists(UniqCaptLR,TotCaptLR,outTab)
% UPDATED 211028 - changed atacHistnorm to HistUniq or HistTot, & changed histogram to histcounts to save time
% Hists for some samples might not get created if samples are almost empty

for iGene = fieldnames(UniqCaptLR)'
    iGene = cell2mat(iGene);
    StartStopArray = [];
    for iSample = fieldnames(UniqCaptLR.(iGene))'
        iSample = cell2mat(iSample);
        StartStopArray=[StartStopArray;
            min(UniqCaptLR.(iGene).(iSample)(:)),max(UniqCaptLR.(iGene).(iSample)(:))]; %finding min/max for each gene and sample
    end
    startPos = min(StartStopArray(:,1)); % min of mins
    stopPos = max(StartStopArray(:,2)); % max of maxs
    StartStop.(iGene)=[startPos,stopPos];
    
    for iSample = fieldnames(UniqCaptLR.(iGene))'
        iSample = cell2mat(iSample);
        
        try % if samples are almost empty, dividing by zero & caculating FragLen won't work
            
            % summing Tn5 insertions per base pair (uniq and tot)
            [UniqBinCounts,~] = histcounts(UniqCaptLR.(iGene).(iSample)(:),startPos:stopPos);
            HistUniqNormLoci.(iGene).(iSample) = UniqBinCounts./(outTab.ACT1_uniq(iSample) + outTab.TUB1_uniq(iSample)...
                + outTab.GAS1_uniq(iSample) + outTab.EFB1_uniq(iSample) + outTab.RPS14B_uniq(iSample));
            HistUniqNormChr.(iGene).(iSample) = UniqBinCounts./(outTab.Chr3_uniq(iSample) + outTab.Chr11_uniq(iSample) + outTab.Chr14_uniq(iSample));
            
            [TotBinCounts,~] = histcounts(TotCaptLR.(iGene).(iSample)(:),startPos:stopPos);
            HistTotNormLoci.(iGene).(iSample) = TotBinCounts./(outTab.ACT1_tot(iSample) + outTab.TUB1_tot(iSample)...
                + outTab.GAS1_tot(iSample) + outTab.EFB1_tot(iSample) + outTab.RPS14B_tot(iSample));
            HistTotNormChr.(iGene).(iSample) = TotBinCounts./(outTab.Chr3_tot(iSample) + outTab.Chr11_tot(iSample) + outTab.Chr14_tot(iSample));
            
            % binning fragment lengths (uniq and tot)
            UniqFragLen = UniqCaptLR.(iGene).(iSample)(:,2) - UniqCaptLR.(iGene).(iSample)(:,1); %Right - Left
            [UniqFragCounts,~] = histcounts(UniqFragLen,[0:10:1000]');
            FragLenHist.uniq.(iGene).(iSample) = UniqFragCounts / length(UniqFragLen); % normalize to probability
            
            TotFragLen = TotCaptLR.(iGene).(iSample)(:,2) - TotCaptLR.(iGene).(iSample)(:,1); %Right - Left
            [TotFragCounts,~] = histcounts(TotFragLen,[0:10:1000]');
            FragLenHist.tot.(iGene).(iSample) = TotFragCounts / length(TotFragLen); % normalize to probability
            
            if ~isempty(UniqFragLen(UniqFragLen>1000))
                warning(['Warning. fragments in ' iGene '_' iSample ' are longer than 1Kb'])
            end
        catch
        end
    end
end