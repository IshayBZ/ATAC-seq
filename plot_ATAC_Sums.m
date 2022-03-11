function SumOut = plot_ATAC_Sums(Exp,datadir,SumIn,AnlOut,Regs)

for iGene = SumIn.Genes
    hfSum.(iGene{1}) = figure('WindowStyle','docked');
    tile.(iGene{1}) = tiledlayout('flow');
    tile.(iGene{1}).TileSpacing = 'compact';
    tile.(iGene{1}).Padding = 'compact';
    for iReg = fieldnames(Regs.(iGene{1}))'
        reg = Regs.(iGene{1}).(iReg{1}) - AnlOut.(iGene{1}).StartStop(1);
        if strcmp(iGene{1},'GAL1')
            MeanON = sum(AnlOut.(iGene{1}).rawHistON.Avg(reg(1):reg(2))) / (reg(2) - reg(1) + 1);
            MeanOFF = sum(AnlOut.(iGene{1}).rawHistOFF.Avg(reg(1):reg(2))) / (reg(2) - reg(1) + 1);
        elseif strcmp(iGene{1},'HO_GAL1')
            MeanON = sum(AnlOut.(iGene{1}).rawHistONSNP.Avg(reg(1):reg(2))) / (reg(2) - reg(1) + 1);
            MeanOFF = sum(AnlOut.(iGene{1}).rawHistOFFSNP.Avg(reg(1):reg(2))) / (reg(2) - reg(1) + 1);
        end
        
        ax = nexttile(tile.(iGene{1}));
        
        SumOut.SumRegTable = table; % collect sums for all conditions (rows) and strains X sorts (columns)
        %     SmplNameCollect = {}; % collect sample names in the same structure of SumRegTable
        for iStrn = SumIn.Strns
            iStrn = {['S_' iStrn{1}]};
            icondsort = 0;
            xlabelCollect.(iStrn{1}) = {};
            for iCond = SumIn.Conds
                for iSort = SumIn.Sorts
                    icondsort = icondsort + 1;
                    CondSortPlotted = 0;
                    try
                        rawHist = AnlOut.(iGene{1}).(iStrn{1}).(iCond{1}).(iSort{1}).rawHist;
                        NAME = nameConvertATACseq_210520(Exp,datadir,[iStrn{1} '_' iCond{1} '_' iSort{1}]);
                        XLABEL = join(NAME(SumIn.Xlabel==1),' _ ');
                        for iBioRep = 1:4
                            if isfield(SumIn,'TechReps')
                                TechReps = SumIn.TechReps;
                            else
                                TechReps = 1:8;
                            end
                            for iTechRep = TechReps
                                if iTechRep==1
                                    MarkFaceCol = [1 1 1];
                                else
                                    MarkFaceCol = SumIn.Colors.(iStrn{1});
                                end
                                try %checking for non-empty reps
                                    SumOut.SumRegTable{iCond{1},[iStrn{1} '_' iSort{1} num2str(iBioRep) '_' num2str(iTechRep)]} ...
                                        = sum(rawHist.Reps{iBioRep,iTechRep}(reg(1):reg(2)));
                                    Mean = sum(rawHist.Reps{iBioRep,iTechRep}(reg(1):reg(2))) / (reg(2) - reg(1) + 1);
                                    if SumIn.Norm==1
                                        plot(icondsort+rand*0.4-0.2,(Mean-MeanOFF)/(MeanON-MeanOFF),'Marker',SumIn.Markers{iBioRep},'MarkerFaceColor',MarkFaceCol,'Color',SumIn.Colors.(iStrn{1}),'LineWidth',2); %SumIn.Colors.(iSort{1})
                                    else
                                        plot(icondsort+rand*0.4-0.2,Mean,'Marker',SumIn.Markers{iBioRep},'MarkerFaceColor',MarkFaceCol,'Color',SumIn.Colors.(iStrn{1}),'LineWidth',2); % SumIn.Colors.(iSort{1})
                                    end
                                    hold on;
                                    CondSortPlotted = 1;
                                catch
                                end
                            end
                        end
                        if CondSortPlotted==1
                            xlabelCollect.(iStrn{1}) = cat(2,xlabelCollect.(iStrn{1}),XLABEL{1});
                        elseif CondSortPlotted==0
                            xlabelCollect.(iStrn{1}) = cat(2,xlabelCollect.(iStrn{1}),' ');
                        end
                    catch
                        xlabelCollect.(iStrn{1}) = cat(2,xlabelCollect.(iStrn{1}),' ');
                    end
                end
            end
        end
        
        %% Configuring Axes and labels
        title([iGene{1} ' - Mean over ' iReg{1}],'Interpreter','none');
        %         xlim([0 icondsort+1]);
        xlim([0 16]);
        if SumIn.Norm==1
            %             ylim([-0.25 2]);
            %             ylim([-0.25 1.5]);
        else
            ylim([0 20E-4]);
        end
        xticks([1:icondsort]);
        xlabelCollect.ALL = xlabelCollect.(iStrn{1}); % uses the last strain in the list
        for iStrn = SumIn.Strns
            iStrn = {['S_' iStrn{1}]};
            for condsort = 1:icondsort
                if ~strcmp(xlabelCollect.(iStrn{1}){condsort},' ')
                    xlabelCollect.ALL{condsort} = xlabelCollect.(iStrn{1}){condsort};
                end
            end
        end
        xticklabels(xlabelCollect.ALL);
        xtickangle(45);
        if strcmp(iReg,'Prom')
            %             legend(SumOut.SumRegTable.Properties.VariableNames,'Interpreter','none','Location','southeast');
        end
    end
    if SumIn.savefigs==1
        saveas(hfSum.(iGene{1}),[SumIn.Dir '/Sum_' iGene{1} '_Norm' num2str(SumIn.Norm) '.png']);
        saveas(hfSum.(iGene{1}),[SumIn.Dir '/Sum_' iGene{1} '_Norm' num2str(SumIn.Norm) '.fig']);
        %             close(hfSum.(iGene{1}));
    end
end
if SumIn.savefigs==1
    CurrentCodeFile = [mfilename '.m'];
    [status,msg] = copyfile(CurrentCodeFile,SumIn.Dir);
    if status==0
        error('Error. Sum-plotting code file was not copied')
    end
    save([SumIn.Dir '/SumIn_' datestr(now,'yymmdd_HHMMSS') '.mat'],'SumIn');
end