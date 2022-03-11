function [rawHist,FragLenHist] = extract_hist(ExtOut,SmplName,iGene,Exp,AnlIn,num2strain,Hist)
% updated 20200130 to save each repeat in its correct column (to track repeat number)
% This function extracts histograms for all repeats of SamplName from HistUniq / HistTot and puts them in a cell array
% and fragment length histograms from ExtOut.FragLenHist.uniq / ExtOut.FragLenHist.tot
% If there is no sample with this name it yields an error message
% updated 211027 to not save FragLen and FragMid into AnlOut - they can be plotted from ExtOut

%% Extract and normalize histograms
%looping thru to find replicate samples
for ibiorep=1:4
    try %looking for SmplName with ibiorep=1 and w/o tech reps
        rawHist.Reps{ibiorep,1} = Hist.(iGene).([SmplName num2str(ibiorep)])';
        rawHist.FullNames{ibiorep,1} = [SmplName num2str(ibiorep)];
        rawHist.SmpNums{ibiorep,1} =  num2strain{strcmp(num2strain{:,2},rawHist.FullNames{ibiorep,1}),1};
        FragLenHist.Reps{ibiorep,1} = ExtOut.FragLenHist.(AnlIn.ut).(iGene).([SmplName num2str(ibiorep)])';
        % Removing samples with low reads
        if ExtOut.outTab{[SmplName num2str(ibiorep)],'SUM_uniq'} < AnlIn.ReadThresh || ismember([SmplName num2str(ibiorep)],AnlIn.BadSmpls)
            rawHist.Bad{ibiorep,1} = rawHist.Reps{ibiorep,1};
            rawHist.Reps{ibiorep,1} = [];
            FragLenHist.Bad{ibiorep,1} = FragLenHist.Reps{ibiorep,1};
            FragLenHist.Reps{ibiorep,1} = [];
            rawHist.AvgTech{ibiorep} = [];
            FragLenHist.AvgTech{ibiorep} = [];
        else
            rawHist.AvgTech{ibiorep} = rawHist.Reps{ibiorep,1};
            FragLenHist.AvgTech{ibiorep} = FragLenHist.Reps{ibiorep,1};
        end
    catch %if smplName w/o tech reps didn't exist look for tech reps
        for itechrep=1:8
            try %
                rawHist.Reps{ibiorep,itechrep} = Hist.(iGene).([SmplName num2str(ibiorep) '_' num2str(itechrep)])';
                rawHist.FullNames{ibiorep,itechrep} = [SmplName num2str(ibiorep) '_' num2str(itechrep)];
                rawHist.SmpNums{ibiorep,itechrep} =  num2strain{strcmp(num2strain{:,2},rawHist.FullNames{ibiorep,itechrep}),1};
                FragLenHist.Reps{ibiorep,itechrep} = ExtOut.FragLenHist.(AnlIn.ut).(iGene).([SmplName num2str(ibiorep) '_' num2str(itechrep)])';
                % Removing samples with low reads
                if ExtOut.outTab{[SmplName num2str(ibiorep) '_' num2str(itechrep)],'SUM_uniq'} < AnlIn.ReadThresh || ismember([SmplName num2str(ibiorep) '_' num2str(itechrep)],AnlIn.BadSmpls)
                    rawHist.Bad{ibiorep,itechrep} = rawHist.Reps{ibiorep,itechrep};
                    rawHist.Reps{ibiorep,itechrep} = [];
                    FragLenHist.Bad{ibiorep,itechrep} = FragLenHist.Reps{ibiorep,itechrep};
                    FragLenHist.Reps{ibiorep,itechrep} = [];
                end
            catch
            end
        end
        try
            rawHist.AvgTech{ibiorep} = mean(cell2mat(rawHist.Reps(ibiorep,:)),2);
            FragLenHist.AvgTech{ibiorep} = mean(cell2mat(FragLenHist.Reps(ibiorep,:)),2);
        catch
        end
    end
end

if exist('rawHist','var')
    rawHist.Name = SmplName;
    rawHist.Avg = mean(cell2mat(rawHist.AvgTech(~cellfun('isempty',rawHist.AvgTech))),2);
    if exist('FragLenHist','var')
        FragLenHist.Avg = mean(cell2mat(FragLenHist.AvgTech(~cellfun('isempty',FragLenHist.AvgTech))),2);
    else
        warning(['Warning. rawHist exists but FragLenHist is missing in condition ' SmplName ' and locus ' iGene]);
        FragLenHist = [];
    end
else
    warning(['Warning. rawHist is missing in condition ' SmplName ' and locus ' iGene]);
    rawHist = [];
    if exist('FragLenHist','var')
        warning(['But FragLenHist exists']);
    else
        warning(['And FragLenHist is missing as well']);
        FragLenHist = [];
    end
end



