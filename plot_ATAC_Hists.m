function plot_ATAC_Hists(Exp,datadir,HistIn,AnlIn,AnlOut,Regs,Sites)
for iGene = HistIn.Genes
    if strcmp(iGene,'HO_GAL1')
        rawHistON = AnlOut.(iGene{1}).rawHistONSNP;
        rawHistOFF = AnlOut.(iGene{1}).rawHistOFFSNP;
    else
        rawHistON = AnlOut.(iGene{1}).rawHistON;
        rawHistOFF = AnlOut.(iGene{1}).rawHistOFF;
    end
    smoothHistON.Avg = smooth(rawHistON.Avg,HistIn.smoothWidth);
    smoothHistOFF.Avg = smooth(rawHistOFF.Avg,HistIn.smoothWidth);
    NAMEon = nameConvertATACseq_210520(Exp,datadir,rawHistON.Name);
    NAMEoff = nameConvertATACseq_210520(Exp,datadir,rawHistOFF.Name);
    LEGENDon = join(NAMEon(HistIn.LegendOnOff==1),' _ ');
    LEGENDoff = join(NAMEoff(HistIn.LegendOnOff==1),' _ ');
    
    StartStopPlot = AnlOut.(iGene{1}).StartStop;
    
    if HistIn.DivByHist==1
        rawHistDiv = AnlOut.(iGene{1}).(HistIn.HistDivStrn).(HistIn.HistDivCond).(HistIn.HistDivSort).rawHist;
        smoothHistDiv.Avg = smooth(rawHistDiv.Avg,HistIn.smoothWidth);
    end
    
    istrncondsort = 0;
    for iStrn = HistIn.Strns
        iStrn = {['S_' iStrn{1}]};
        labelCollect={}; % collect labels for all conditions of iStrn
        hfHist.(iStrn{1}) = figure('WindowStyle','docked');
        tile.(iStrn{1}) = tiledlayout('flow');
        tile.(iStrn{1}).TileSpacing = 'compact';
        tile.(iStrn{1}).Padding = 'compact';
        ax = nexttile(tile.(iStrn{1}));
        
        % PLOT mean ON,OFF HISTOGRAMS
        if HistIn.plotOnOffHists==1
            labelCollect = cat(2,labelCollect,LEGENDon{1});
            labelCollect = cat(2,labelCollect,LEGENDoff{1});
            if HistIn.DivByHist==0
                plot(StartStopPlot(1):StartStopPlot(2)-1,smoothHistOFF.Avg,'--','LineWidth',1.5,'Color',[0 0 1]); %[0.93 0.69 0.13]
                hold on
                plot(StartStopPlot(1):StartStopPlot(2)-1,smoothHistON.Avg,'-','LineWidth',1.5,'Color',[0 0 1]);
            elseif HistIn.DivByHist==1
                plot(StartStopPlot(1):StartStopPlot(2)-1,smoothHistOFF.Avg./smoothHistDiv.Avg,'--','LineWidth',1.5,'Color',[0 0 1]);
                hold on
                plot(StartStopPlot(1):StartStopPlot(2)-1,smoothHistON.Avg./smoothHistDiv.Avg,'-','LineWidth',1.5,'Color',[0 0 1]);
            end
        end
        icond = 0;
        for iCond = HistIn.Conds
            if isfield(HistIn,'Sorts')
                Sorts = HistIn.Sorts;
            else
                Sorts = fieldnames(AnlOut.(iGene{1}).(iStrn{1}).(iCond{1}))';
            end
            isort = 0;
            icond = icond + 1;
            for iSort = Sorts
                istrncondsort = istrncondsort + 1;
                isort = isort + 1;
                if HistIn.ColorBySort==1
                    icolor = isort;
                    style = HistIn.CondStyles{icond};
                elseif HistIn.ColorByCond==1
                    icolor = icond;
                    style = HistIn.SortStyles{isort};
                else
                    icolor = istrncondsort;
                end
                try
                    rawHist = AnlOut.(iGene{1}).(iStrn{1}).(iCond{1}).(iSort{1}).rawHist;
                    smoothHist.Avg = smooth(rawHist.Avg,HistIn.smoothWidth);
%                     [STRN,COND,TAGMENT,CELLNUM] = nameConvertATACseq_210520(Exp,datadir,rawHist.Name);
                    NAME = nameConvertATACseq_210520(Exp,datadir,rawHist.Name);
                    LEGEND = join(NAME(HistIn.Legend==1),' _ ');
                    
                    if HistIn.AvgRep==1
                        labelCollect = cat(2,labelCollect,LEGEND{1});
                        if HistIn.DivByHist==0
                            plot(StartStopPlot(1):StartStopPlot(2)-1,smoothHist.Avg,'LineWidth',1.5,'Color',HistIn.Colors{icolor},'LineStyle',style);
                        else
                            plot(StartStopPlot(1):StartStopPlot(2)-1,smoothHist.Avg./smoothHistDiv.Avg,'LineWidth',1.5,'Color',HistIn.Colors{icolor},'LineStyle',style);
                        end
                        hold on
                    elseif HistIn.AvgRep==0 && HistIn.AvgTechRep==1
                        BioReps = find(~cellfun(@isempty,AnlOut.(iGene{1}).(iStrn{1}).(iCond{1}).(iSort{1}).rawHist.AvgTech));
                        for ibiorep = BioReps
                            smoothHist.BioReps{ibiorep} = smooth(rawHist.AvgTech{ibiorep},HistIn.smoothWidth);
                            labelCollect = cat(2,labelCollect,[LEGEND{1} ' _ ' num2str(ibiorep)]);
                            if ~isempty(smoothHist.BioReps{ibiorep}) % because some reps might be empty
                                if HistIn.DivByHist==0
                                    plot(StartStopPlot(1):StartStopPlot(2)-1,smoothHist.BioReps{ibiorep},'LineWidth',1.5,'Color',HistIn.Colors{icolor},'LineStyle',style);
                                else
                                    plot(StartStopPlot(1):StartStopPlot(2)-1,smoothHist.BioReps{ibiorep}./smoothHistDiv.Avg,'LineWidth',1.5,'Color',HistIn.Colors{icolor},'LineStyle',style);
                                end
                                hold on;
                            end
                        end
                    elseif HistIn.AvgRep==0 && HistIn.AvgTechRep==0
                        if isfield(HistIn,'Reps')
                            Reps = HistIn.Reps;
                        else
                            Reps = find(~cellfun(@isempty,AnlOut.(iGene{1}).(iStrn{1}).(iCond{1}).(iSort{1}).rawHist.Reps));
                            Reps = reshape(Reps,[1,length(Reps)]); % Turn Reps (which could be a row/column) into a row
                        end
                        for irep = Reps % linear indices of bio X tech repeats
                            [ibiorep,itechrep] = ind2sub(size(AnlOut.(iGene{1}).(iStrn{1}).(iCond{1}).(iSort{1}).rawHist.Reps),irep);
                            smoothHist.Reps{ibiorep,itechrep} = smooth(rawHist.Reps{ibiorep,itechrep},HistIn.smoothWidth);
                            labelCollect = cat(2,labelCollect,[LEGEND{1} ' _ ' num2str(ibiorep) '_' num2str(itechrep)]);

                            if ~isempty(smoothHist.Reps{ibiorep,itechrep}) % because some reps might be empty
                                if HistIn.DivByHist==0
                                    plot(StartStopPlot(1):StartStopPlot(2)-1,smoothHist.Reps{ibiorep,itechrep},'LineWidth',1.5,'Color',HistIn.Colors{icolor},'LineStyle',style);
                                else
                                    plot(StartStopPlot(1):StartStopPlot(2)-1,smoothHist.Reps{ibiorep,itechrep}./smoothHistDiv.Avg,'LineWidth',1.5,'Color',HistIn.Colors{icolor},'LineStyle',style);
                                end
                                hold on;
                            end
                        end
                    end
                catch
                end
            end
        end
        
        %% Drawing genomic features
        if HistIn.DivByHist==0
            if strncmp(iGene{1},'ZF',2)
                rectHeight=0.02e-3;
                rectY=-0.044e-3;
            else
                rectHeight=0.01e-3;
                rectY=-0.022e-3;
            end
        else
            rectHeight=0.01;
            rectY=-0.022;
        end
        CDS(1) = min(AnlOut.CoordsTable{iGene{1},'CDS_5'},AnlOut.CoordsTable{iGene{1},'CDS_3'});
        CDS(2) = max(AnlOut.CoordsTable{iGene{1},'CDS_5'},AnlOut.CoordsTable{iGene{1},'CDS_3'});
        rectangle('Position',[max(CDS(1),StartStopPlot(1)),rectY,min(CDS(2),StartStopPlot(2)) - max(CDS(1),StartStopPlot(1)),rectHeight],'FaceColor',[0 0 0]);
        annotation('textbox',[0.15 0.75 0.2 0.05],'String',strcat(AnlOut.CoordsTable{iGene{1},'CDS'},' CDS'),'Color',[0 0 0],'EdgeColor','none','FontSize',14);

        if ~isnan(AnlOut.CoordsTable{iGene{1},'U1_CDS_5'}) && ~isempty(AnlOut.CoordsTable{iGene{1},'U1_CDS_5'})
            U1_CDS(1) = min(AnlOut.CoordsTable{iGene{1},'U1_CDS_5'},AnlOut.CoordsTable{iGene{1},'U1_CDS_3'});
            U1_CDS(2) = max(AnlOut.CoordsTable{iGene{1},'U1_CDS_5'},AnlOut.CoordsTable{iGene{1},'U1_CDS_3'});
            try
                rectangle('Position',[max(U1_CDS(1),StartStopPlot(1)),rectY+rectHeight/2,min(U1_CDS(2),StartStopPlot(2)) - max(U1_CDS(1),StartStopPlot(1)),rectHeight],'FaceColor',[0 0 1]);
                annotation('textbox',[0.15 0.8 0.2 0.05],'String',strcat(AnlOut.CoordsTable{iGene{1},'U1_CDS'},' CDS'),'Color',[0 0 1],'EdgeColor','none','FontSize',14);
            catch end
        end
        
        %% Highlight summing regions
        if isfield(Regs,iGene{1})
            for iReg=fieldnames(Regs.(iGene{1}))'
                if strcmp(iReg{1},'Prom')
                    continue
                end
                if strcmp(iReg{1},'NFR')
                    rectangle('Position',[Regs.(iGene{1}).(iReg{1})(1),rectY-rectHeight,Regs.(iGene{1}).(iReg{1})(2) - Regs.(iGene{1}).(iReg{1})(1),rectHeight],'FaceColor',[0.7 0.7 0.7]);
                else
                    rectangle('Position',[Regs.(iGene{1}).(iReg{1})(1),rectY-rectHeight,Regs.(iGene{1}).(iReg{1})(2) - Regs.(iGene{1}).(iReg{1})(1),rectHeight],'FaceColor',[0.4 0.4 0.4]);
                end
            end
        end
        
        %% Highlight promoter features - FIX LATER TO ACCOUNT FOR FEATURES DELETED FROM PROMOTER
        try
            tata(1) = min(Sites.TATA.(iGene{1})(1),Sites.TATA.(iGene{1})(2));
            tata(2) = max(Sites.TATA.(iGene{1})(1),Sites.TATA.(iGene{1})(2));
            rectangle('Position',[tata(1),rectY,tata(2)-tata(1)+1,rectHeight],'FaceColor',[0 1 0]);
        catch end
        try
            tss(1) = min(Sites.TSS.(iGene{1})(1),Sites.TSS.(iGene{1})(2));
            tss(2) = max(Sites.TSS.(iGene{1})(1),Sites.TSS.(iGene{1})(2));
            rectangle('Position',[tss(1),rectY,tss(2)-tss(1)+1,rectHeight],'FaceColor',[1 0 1]);
        catch end
        try
            for iGal4 = 1:4
                gal4bs(1) = min(Sites.Gal4bs.(iGene{1}){iGal4}(1),Sites.Gal4bs.(iGene{1}){iGal4}(2));
                gal4bs(2) = max(Sites.Gal4bs.(iGene{1}){iGal4}(1),Sites.Gal4bs.(iGene{1}){iGal4}(2));
                rectangle('Position',[gal4bs(1),rectY,gal4bs(2)-gal4bs(1)+1,rectHeight],'FaceColor',[0 0 1]);
            end
        catch end
        try
            for iMig1 = 1:2
                mig1bs(1) = min(Sites.Mig1bs.(iGene{1}){iMig1}(1),Sites.Mig1bs.(iGene{1}){iMig1}(2));
                mig1bs(2) = max(Sites.Mig1bs.(iGene{1}){iMig1}(1),Sites.Mig1bs.(iGene{1}){iMig1}(2));
                rectangle('Position',[mig1bs(1),rectY,mig1bs(2)-mig1bs(1)+1,rectHeight],'FaceColor',[1 0 0]);
            end
        catch end
        try
            core(1) = min(Sites.CORE.(iGene{1})(1),Sites.CORE.(iGene{1})(2));
            core(2) = max(Sites.CORE.(iGene{1})(1),Sites.CORE.(iGene{1})(2));
            rectangle('Position',[core(1),rectY,core(2)-core(1)+1,rectHeight*3],'FaceColor',[0 1 1]);
        catch end
        try
            ter(1) = min(Sites.TER.(iGene{1})(1),Sites.TER.(iGene{1})(2));
            ter(2) = max(Sites.TER.(iGene{1})(1),Sites.TER.(iGene{1})(2));
            rectangle('Position',[ter(1),rectY,ter(2)-ter(1)+1,rectHeight*3],'FaceColor',[1 1 0]);
        catch end
        try
            for iZF = 1:2
                zfbs(1) = min(Sites.ZFbs.(iGene{1}){iZF}(1),Sites.ZFbs.(iGene{1}){iZF}(2));
                zfbs(2) = max(Sites.ZFbs.(iGene{1}){iZF}(1),Sites.ZFbs.(iGene{1}){iZF}(2));
                rectangle('Position',[zfbs(1),rectY,zfbs(2)-zfbs(1)+1,rectHeight*2],'FaceColor',[0 0 1]);
            end
        catch end
        
        %% Configuring Axes and labels
        title([Exp ' _ ' iGene{1} ' _ ' NAME{2}],'Interpreter','none');
        if HistIn.plotlegend==1
            legend(labelCollect,'FontSize',10,'Location','eastoutside','Interpreter','none');
        end
        
        if U1_CDS(2)<CDS(1)
            xlim([max(U1_CDS(2)-200,StartStopPlot(1)), min(CDS(1)+200,StartStopPlot(2))]);
        elseif CDS(2)<U1_CDS(1)
            xlim([max(CDS(2)-200,StartStopPlot(1)), min(U1_CDS(1)+200,StartStopPlot(2))]);
            if strcmp(iGene{1},'HO_GAL1')
%                 rectangle('Position',[CDS(2)+158,rectY-rectHeight/2,8,rectHeight*3],'FaceColor',[0.5 0 0]);
                xlim([max(CDS(2)-200,StartStopPlot(1)), StartStopPlot(2)-500]); %StartStopPlot(2)-125
                set(gca,'XDir','reverse');
            end
        end
        if HistIn.DivByHist == 0
            if strcmp(AnlIn.ut,'uniq')
                if strncmp(iStrn{1},'S_ZF',4)
                    ylim([-7E-5 26E-4]);
                else
                    ylim([-5E-5 12E-4]);
                end
            elseif strcmp(AnlIn.ut,'tot')
                ylim([-5E-5 2E-3]);
            end
        elseif HistIn.DivByHist == 1
            
        end
        
        if HistIn.savefigs==1
            saveas(hfHist.(iStrn{1}),[HistIn.Dir '/Hist_' iGene{1} '_' iStrn{1} '_Div' '.png']);
            saveas(hfHist.(iStrn{1}),[HistIn.Dir '/Hist_' iGene{1} '_' iStrn{1} '_Div' '.fig']);
        end
    end
end
if HistIn.savefigs==1
    CurrentCodeFile = [mfilename '.m'];
    [status,msg] = copyfile(CurrentCodeFile,HistIn.Dir);
    if status==0
        error('Error. hist-plotting code file was not copied')
    end
    save([HistIn.Dir '/HistIn_' datestr(now,'yymmdd_HHMMSS') '.mat'],'HistIn');
end