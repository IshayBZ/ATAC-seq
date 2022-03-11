% updated 20200130 to arrange only (no plotting or summing) AND to track repeat number
% GOAL: arrange all samples of a given eperiment in a strcut saved in an Analysis folder together with a code snapshot

function [AnlDir,AnlOut] = a_Analyze_ATAC(Exp,datadir,ExtDir,ExtOut,Hist,AnlIn,SmplNameON,SmplNameOFF,SmplNameONSNP,SmplNameOFFSNP,Regs)
%% Read CDS coordinates from IDglobal/Coords sheet
try
    AnlOut.CoordsTable = readtable([datadir '../../IDglobal_' Exp '.xlsx'],'ReadRowNames',true,'Sheet', 'Coords');
catch
    AnlOut.CoordsTable = readtable([datadir 'IDglobal_' Exp '.xlsx'],'ReadRowNames',true,'Sheet', 'Coords');
end

%% Read original sample numbers from IDglobal/Samples sheet
try
    AnlOut.num2strain = readtable([datadir '../../IDglobal_',Exp,'.xlsx'],'Sheet','Samples');
catch
    AnlOut.num2strain = readtable([datadir 'IDglobal_',Exp,'.xlsx'],'Sheet','Samples');
end
NumSmpls = size(AnlOut.num2strain,1);

%% Arrange samples in out struct, based on the 'Table' sheet in IDglobal_L93.xlsx
try
    AnlOut.SmplTable = readtable([datadir '../../IDglobal_' Exp '.xlsx'],'ReadRowNames',true,'Sheet', 'Table');
catch
    AnlOut.SmplTable = readtable([datadir 'IDglobal_' Exp '.xlsx'],'ReadRowNames',true,'Sheet', 'Table');
end

% for iGene = fieldnames(ExtOut.StartStop)'
for iGene = fieldnames(Hist)'
    AnlOut.(iGene{1}).StartStop = ExtOut.StartStop.(iGene{1}); % save start-stop positions to Analysis struct
    
    %% Extracting scaling histograms (OFF and ON samples)
    if strcmp(iGene{1},'HO_GAL1') % 220207 change to extract rawhistONSNP only for HO_GAL1, and don't extract it for other genes
        [AnlOut.(iGene{1}).rawHistONSNP,~] = extract_hist_211027(ExtOut,SmplNameONSNP,iGene{1},Exp,AnlIn,AnlOut.num2strain,Hist);
        [AnlOut.(iGene{1}).rawHistOFFSNP,~] = extract_hist_211027(ExtOut,SmplNameOFFSNP,iGene{1},Exp,AnlIn,AnlOut.num2strain,Hist);
    else
        [AnlOut.(iGene{1}).rawHistON,~] = extract_hist_211027(ExtOut,SmplNameON,iGene{1},Exp,AnlIn,AnlOut.num2strain,Hist);
        [AnlOut.(iGene{1}).rawHistOFF,~] = extract_hist_211027(ExtOut,SmplNameOFF,iGene{1},Exp,AnlIn,AnlOut.num2strain,Hist);
    end
    
    %% Extracting all histograms (for all strain-condition combinations tested)
    AnlOut.SmplCount = 0; % counting samples to make sure we got them all
    StrnNum = length(AnlOut.SmplTable.Properties.VariableNames);
    for iStrn = AnlOut.SmplTable.Properties.VariableNames
        labelCollect={}; % collect labels for all samples of a given strain
        
        for iCond = AnlOut.SmplTable.Properties.RowNames' % doesn't work w/o transpose that makes it a 'row-vector' cell array
            if ~isempty(AnlOut.SmplTable{iCond{1},iStrn{1}}{1}) % this char vector notes all sort gates (other params in general, e.g. Tagmentation method) for this strain-cond combination separated by ',' (e.g. 'F,O')
                Sorts = AnlOut.SmplTable{iCond{1},iStrn{1}}{1};
                Sorts = regexp(Sorts,',','split');
                for iSort = Sorts %looping thru sort gates (params in general)
                    SmplName = [iStrn{1} '_' iCond{1} '_' iSort{1}]; % this name doesn't include repeat number and Exp number
                    
                    % Fetch Histograms for all repeats of SmplName - i.e. samples whose name start with SmplName
                    [rawHist,FragLenHist] = extract_hist_211027(ExtOut,SmplName,iGene{1},Exp,AnlIn,AnlOut.num2strain,Hist);
                    % rawHist.(iSort{1}) is a struct with Name,Avg,Reps,Bad fields, Reps and Bad are cell arrays of good and bad reps (inc empty cells for nonexisting reps)
                    
                    Reps = find(~cellfun(@isempty,rawHist.Reps)); % finds actual (non-empty) reps
                    AnlOut.SmplCount = AnlOut.SmplCount + length(Reps);
                    
                    AnlOut.(iGene{1}).(iStrn{1}).(iCond{1}).(iSort{1}).rawHist = rawHist;
                    AnlOut.(iGene{1}).(iStrn{1}).(iCond{1}).(iSort{1}).FragLenHist = FragLenHist;
                    
                    %% extract accessibility distance from average for qc plot
                    % use either ZF or GAL1 gene
                    if strncmp(iStrn{1},'S_ZF1',5)
                        if strcmp(iGene{1},'ZF1_CYC1')
                            reg = Regs.(iGene{1}).Nuc_m1 - AnlOut.(iGene{1}).StartStop(1);
                            Reps = reshape(Reps,[1,length(Reps)]); % Turn Reps (which could be a row/column) into a row
                            for irep = Reps
                                [ibiorep,itechrep] = ind2sub(size(rawHist.Reps),irep);
                                Acc = sum(rawHist.Reps{ibiorep,itechrep}(reg(1):reg(2))) / (reg(2) - reg(1) + 1);
                                AccAvg = sum(rawHist.Avg(reg(1):reg(2))) / (reg(2) - reg(1) + 1);
                                RepName = rawHist.FullNames{ibiorep,itechrep};
                                AnlOut.Dist.(RepName) = abs(Acc-AccAvg)/AccAvg;
                            end
                        end
                    elseif strncmp(iStrn{1},'S_ZF2',5)
                        if strcmp(iGene{1},'ZF2_CYC1')
                            reg = Regs.(iGene{1}).Nuc_m1 - AnlOut.(iGene{1}).StartStop(1);
                            Reps = reshape(Reps,[1,length(Reps)]); % Turn Reps (which could be a row/column) into a row
                            for irep = Reps
                                [ibiorep,itechrep] = ind2sub(size(rawHist.Reps),irep);
                                Acc = sum(rawHist.Reps{ibiorep,itechrep}(reg(1):reg(2))) / (reg(2) - reg(1) + 1);
                                AccAvg = sum(rawHist.Avg(reg(1):reg(2))) / (reg(2) - reg(1) + 1);
                                RepName = rawHist.FullNames{ibiorep,itechrep};
                                AnlOut.Dist.(RepName) = abs(Acc-AccAvg)/AccAvg;
                            end
                        end
                    else %for all other strains, use GAL1
                        if strcmp(iGene{1},'GAL1')
                            reg = Regs.(iGene{1}).Nuc_m1 - AnlOut.(iGene{1}).StartStop(1);
                            Reps = reshape(Reps,[1,length(Reps)]); % Turn Reps (which could be a row/column) into a row
                            for irep = Reps
                                [ibiorep,itechrep] = ind2sub(size(rawHist.Reps),irep);
                                Acc = sum(rawHist.Reps{ibiorep,itechrep}(reg(1):reg(2))) / (reg(2) - reg(1) + 1);
                                AccAvg = sum(rawHist.Avg(reg(1):reg(2))) / (reg(2) - reg(1) + 1);
                                RepName = rawHist.FullNames{ibiorep,itechrep};
                                AnlOut.Dist.(RepName) = abs(Acc-AccAvg)/AccAvg;
                            end
                        end
                    end
                end
            end
        end
    end
end
if AnlIn.checkSmplCount==1 && ExtOut.SmplCount ~= AnlOut.SmplCount
    %     error('Error. Table and Samples sheets do not match')
end

%% Calculate & plot read statistics
% icds = 0;
AnlOut.CaptLength = 0;
for iCDS = {'ACT1','CAT8','EFB1','HO_GAL1','GAL1','GAL2','GAL3','GAL4','GAL7','GAL80','GAS1','GCY1','HXT1','MAL31','MAL33','RPS14B','SUC2','TUB1'}
    %     icds = icds +1;
    try
        AnlOut.CaptLength = AnlOut.CaptLength + AnlOut.CoordsTable{iCDS,'Capt_stop'} - AnlOut.CoordsTable{iCDS,'Capt_start'};
        %     AnlOut.CaptLength = AnlOut.CaptLength + AnlOut.CoordsTable.Capt_stop(icds) - AnlOut.CoordsTable.Capt_start(icds);
    catch
        warning([iCDS{1} ' is missing in the Coords tab of IDglobal and was dropped from the capture length calculation']);
    end
    % overall length of captured loci
end

hf = figure('WindowStyle','docked');
tile = tiledlayout(3,1); tile.TileSpacing = 'compact'; tile.Padding = 'compact';

% calculate averages
%Coverage
AnlOut.uniqCAPTavg = sum(ExtOut.outTab.SUM_uniq)/NumSmpls;
AnlOut.totCAPTavg = sum(ExtOut.outTab.SUM_tot)/NumSmpls;
AnlOut.uniqALLavg = sum(ExtOut.outTab.SUMCHR_uniq)/NumSmpls;
AnlOut.totALLavg = sum(ExtOut.outTab.SUMCHR_tot)/NumSmpls;

%Duplication
AnlOut.dupCAPTavg = AnlOut.totCAPTavg / AnlOut.uniqCAPTavg;
AnlOut.dupUNCAPTavg = (AnlOut.totALLavg - AnlOut.totCAPTavg) / (AnlOut.uniqALLavg - AnlOut.uniqCAPTavg);

% Enrichment
AnlOut.totCAPTperBPavg = AnlOut.totCAPTavg / AnlOut.CaptLength;
AnlOut.totUNCAPTperBPavg = (AnlOut.totALLavg - AnlOut.totCAPTavg) / (12000000 - AnlOut.CaptLength);
AnlOut.CaptENRICHavg = AnlOut.totCAPTperBPavg / AnlOut.totUNCAPTperBPavg;

for iSmpl = ExtOut.outTab.Properties.RowNames'
    Smplnum = AnlOut.num2strain{strcmp(AnlOut.num2strain{:,2},iSmpl{1}),1};
    
    %%% Calculate coverage
    AnlOut.uniqCAPT.(iSmpl{1}) = ExtOut.outTab.SUM_uniq(iSmpl{1});
    AnlOut.totCAPT.(iSmpl{1}) = ExtOut.outTab.SUM_tot(iSmpl{1});
    AnlOut.uniqALL.(iSmpl{1}) = ExtOut.outTab.SUMCHR_uniq(iSmpl{1});
    AnlOut.totALL.(iSmpl{1}) = ExtOut.outTab.SUMCHR_tot(iSmpl{1});
    % Plot
    ax1 = nexttile(1);
    plot(Smplnum,AnlOut.uniqCAPT.(iSmpl{1}),'bo',Smplnum,AnlOut.totCAPT.(iSmpl{1}),'ko'); %,Smplnum,AnlOut.uniqALL.(iSmpl{1}),'bd',Smplnum,AnlOut.totALL.(iSmpl{1}),'kd');
    hold on;
    if AnlOut.uniqALL.(iSmpl{1}) < AnlIn.ReadThresh
        plot(Smplnum,AnlOut.uniqCAPT.(iSmpl{1}),'rx',Smplnum,AnlOut.totCAPT.(iSmpl{1}),'rx'); %,Smplnum,AnlOut.uniqALL.(iSmpl{1}),'rx',Smplnum,AnlOut.totALL.(iSmpl{1}),'rx');
    end
    
    %%% Calculate duplication
    AnlOut.dupCAPT.(iSmpl{1}) = ExtOut.outTab.SUM_tot(iSmpl{1})/ExtOut.outTab.SUM_uniq(iSmpl{1});
    AnlOut.dupUNCAPT.(iSmpl{1}) = (AnlOut.totALL.(iSmpl{1}) - AnlOut.totCAPT.(iSmpl{1})) / (AnlOut.uniqALL.(iSmpl{1}) - AnlOut.uniqCAPT.(iSmpl{1}));
    % Plot
    ax2 = nexttile(2);
    plot(Smplnum,AnlOut.dupCAPT.(iSmpl{1}),'mo',Smplnum,AnlOut.dupUNCAPT.(iSmpl{1}),'co');
    hold on;
    if AnlOut.uniqCAPT.(iSmpl{1}) < AnlIn.ReadThresh
        plot(Smplnum,AnlOut.dupCAPT.(iSmpl{1}),'rx',Smplnum,AnlOut.dupUNCAPT.(iSmpl{1}),'rx');
    end
    
    % Calculate enrichment
    AnlOut.totCAPTperBP.(iSmpl{1}) = AnlOut.totCAPT.(iSmpl{1}) / AnlOut.CaptLength;
    AnlOut.totUNCAPTperBP.(iSmpl{1}) = (AnlOut.totALL.(iSmpl{1}) - AnlOut.totCAPT.(iSmpl{1})) / (12000000 - AnlOut.CaptLength);
    AnlOut.CaptENRICH.(iSmpl{1}) = AnlOut.totCAPTperBP.(iSmpl{1}) / AnlOut.totUNCAPTperBP.(iSmpl{1});
    
    ax3 = nexttile(3);
    plot(Smplnum,AnlOut.CaptENRICH.(iSmpl{1}),'ko');
    hold on;
    title(['Enrichment = captured read per bp / uncaptured read per bp,       Avg enrichment = ' num2str(AnlOut.CaptENRICHavg,3)]);
    
    %     legend(captured in black, ALL in red);
end
xlabel(tile,'Sample number');

axes(ax1); ax1.YScale = 'log'; ylabel('Coverage'); ylim([AnlIn.ReadThresh 30000000]);
title(['Coverage = # of reads,     Avg capt reads (uniq,b) = ' num2str(AnlOut.uniqCAPTavg,2)...
    '    Avg capt reads (tot,k) = ' num2str(AnlOut.totCAPTavg,2) '     per bp (depth) = ' num2str(AnlOut.totCAPTperBPavg,2) ]);

axes(ax2); ax2.YScale = 'log'; ylabel('PCR duplication');
title(['PCR Duplication = Tot reads / Uniq reads,     Avg Dup (capt,m) = ' num2str(AnlOut.dupCAPTavg,2) '     Norm by per sample depth = ' num2str(AnlOut.dupCAPTavg/AnlOut.totCAPTperBPavg,2)...
    '   Avg Dup (uncapt,c) = '   num2str(AnlOut.dupUNCAPTavg,2) '     Norm by per sample depth = ' num2str(AnlOut.dupUNCAPTavg/AnlOut.totUNCAPTperBPavg,2)]);

axes(ax3); ylabel('Enrichment'); ylim([0 600]);

%% plot qc - duplication / read count vs noise (distance from average)
hf2 = figure('WindowStyle','docked');
tile = tiledlayout(2,2); tile.TileSpacing = 'compact'; tile.Padding = 'compact';

DIST = []; DUPCAPT = []; UNIQCAPT = []; DNACONC = [];
for iSmpl = fieldnames(AnlOut.Dist)' % goes over all samples except those in BAD
    DIST = [DIST AnlOut.Dist.(iSmpl{1})];
    DUPCAPT = [DUPCAPT AnlOut.dupCAPT.(iSmpl{1})];
    UNIQCAPT = [UNIQCAPT AnlOut.uniqCAPT.(iSmpl{1})];
    AnlOut.DNAconc.(iSmpl{1}) = AnlOut.num2strain{strcmp(AnlOut.num2strain{:,2},iSmpl{1}),3};
    DNACONC = [DNACONC AnlOut.DNAconc.(iSmpl{1})];
end

ax1 = nexttile(1);
plot(DIST,DUPCAPT,'ko');
[R1,P1] = corrcoef(DIST,DUPCAPT);
axes(ax1); ax1.YScale = 'log'; ylabel('PCR duplication'); xlabel('Distance from average');
title(['R = ' num2str(R1(1,2),2) '   P = ' num2str(P1(1,2),3)]);

ax2 = nexttile(2);
plot(DIST,UNIQCAPT,'ko');
[R2,P2] = corrcoef(DIST,UNIQCAPT);
axes(ax2); ax2.YScale = 'log'; ylabel('# of unique reads'); xlabel('Distance from average');
title(['R = ' num2str(R2(1,2),2) '   P = ' num2str(P2(1,2),3)]);

% PLOT OF DNA CONC
ax3 = nexttile(3);
plot(DIST,DNACONC,'ko');
[R3,P3] = corrcoef(DIST,DNACONC);
axes(ax3); ylabel('DNA concentration [ug/ml]'); xlabel('Distance from average');
title(['R = ' num2str(R3(1,2),2) '   P = ' num2str(P3(1,2),3)]);

ax4 = nexttile(4);
plot(DUPCAPT,DNACONC,'ko');
[R4,P4] = corrcoef(DUPCAPT,DNACONC);
axes(ax4); ylabel('DNA concentration [ug/ml]'); xlabel('PCR duplication');
title(['R = ' num2str(R4(1,2),2) '   P = ' num2str(P4(1,2),3)]);

%% save data, stat fig, and code snapshot
AnlDir = [ExtDir '/Anl_' datestr(now,'yymmdd_HHMM') '_' AnlIn.ut '_' AnlIn.Norm '_Th' num2str(AnlIn.ReadThresh/1000) 'k']; % to save the analysis mat file
mkdir(AnlDir);
CurrentCodeFile = [mfilename '.m'];
[status,msg] = copyfile(CurrentCodeFile,AnlDir);
if status==0
    error('Error. File was not copied')
end
save([AnlDir '/AnlOut.mat'],'AnlOut');
saveas(hf,[AnlDir '/Summary_stats.png']);
saveas(hf,[AnlDir '/Summary_stats.fig']);
saveas(hf2,[AnlDir '/QC.png']);
saveas(hf2,[AnlDir '/QC.fig']);