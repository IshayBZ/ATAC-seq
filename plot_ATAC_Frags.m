function plot_ATAC_Frags(Exp,FragIn,AnlOut,Regs,Sites)
for iGene = FragIn.Genes
    StartStopPlot = AnlOut.(iGene{1}).StartStop;
    
    if FragIn.stackStrns==1
        labelCollect={}; % collect labels for all conditions of iStrn
        hfLenHist.(iGene{1}) = figure('WindowStyle','docked');
        tile.(iGene{1}) = tiledlayout('flow');
        tile.(iGene{1}).TileSpacing = 'compact';
        tile.(iGene{1}).Padding = 'compact';
        ax = nexttile(tile.(iGene{1}));
    end
    istrncondsort = 0;
    for iStrn = FragIn.Strns
        iStrn = {['S_' iStrn{1}]};
        if FragIn.stackStrns==0
            labelCollect={}; % collect labels for all conditions of iStrn
            hfLenHist.(iGene{1}).(iStrn{1}) = figure('WindowStyle','docked'); %,'Position',[]
            tile.(iGene{1}).(iStrn{1}) = tiledlayout('flow');
            tile.(iGene{1}).(iStrn{1}).TileSpacing = 'compact';
            tile.(iGene{1}).(iStrn{1}).Padding = 'compact';
            ax = nexttile(tile.(iGene{1}).(iStrn{1}));
        end
        
        icond = 0;
        for iCond = FragIn.Conds
            if isfield(FragIn,'Sorts')
                Sorts = FragIn.Sorts;
            else
                Sorts = fieldnames(AnlOut.(iGene{1}).(iStrn{1}).(iCond{1}))';
            end
            isort = 0;
            icond = icond + 1;
            for iSort = Sorts
                istrncondsort = istrncondsort + 1;
                isort = isort + 1;
                if FragIn.ColorBySort==1
                    icolor = isort;
                else
                    icolor = istrncondsort;
                end
                try
                    FragLenHist = AnlOut.(iGene{1}).(iStrn{1}).(iCond{1}).(iSort{1}).FragLenHist;
                    rawHist = AnlOut.(iGene{1}).(iStrn{1}).(iCond{1}).(iSort{1}).rawHist;
                    FragLen = AnlOut.(iGene{1}).(iStrn{1}).(iCond{1}).(iSort{1}).FragLen;
                    FragMid = AnlOut.(iGene{1}).(iStrn{1}).(iCond{1}).(iSort{1}).FragMid;
                    
                    if FragIn.plotlenVSmid == 1 %plotting fragment size vs midpoint of chosen rep of gluc sample
                        if icond>1
                            ax = nexttile(tile.(iGene{1}).(iStrn{1})); % to not overlay different conds on the same V plot
                        end
                        labelCollect = cat(2,labelCollect,[Exp '_' rawHist.Name]);
                        plot(FragMid.Reps{FragIn.lenVSmidBioRep,FragIn.lenVSmidTechRep},FragLen.Reps{FragIn.lenVSmidBioRep,FragIn.lenVSmidTechRep},'o','MarkerSize',1);
                        ylim([0 400]);
                    else % plotting Fragment length distribution
                        if FragIn.AvgRep==1
                            labelCollect = cat(2,labelCollect,[Exp '_' rawHist.Name]);
                            plot([10:10:1000],FragLenHist.Avg,'LineWidth',1.5,'Color',FragIn.Colors{icolor},'LineStyle',FragIn.CondStyles{icond});
                            hold on
                        elseif FragIn.AvgRep==0 && FragIn.AvgTechRep==1
                            BioReps = find(~cellfun(@isempty,AnlOut.(iGene{1}).(iStrn{1}).(iCond{1}).(iSort{1}).FragLenHist.AvgTech));
                            for ibiorep = BioReps
                                labelCollect = cat(2,labelCollect,[Exp '_' rawHist.Name num2str(ibiorep)]);
                                if ~isempty(FragLenHist.BioReps{ibiorep}) % because some reps might be empty
                                    plot([10:10:1000],FragLenHist.BioReps{ibiorep},'LineWidth',1.5,'Color',FragIn.Colors{icolor},'LineStyle',FragIn.CondStyles{icond});
                                    hold on;
                                end
                            end
                        elseif FragIn.AvgRep==0 && FragIn.AvgTechRep==0
                            if isfield(FragIn,'Reps')
                                Reps = FragIn.Reps;
                            else
                                Reps = find(~cellfun(@isempty,AnlOut.(iGene{1}).(iStrn{1}).(iCond{1}).(iSort{1}).FragLenHist.Reps));
                                Reps = reshape(Reps,[1,length(Reps)]); % Turn Reps (which could be a row/column) into a row
                            end
                            for irep = Reps % linear indices of bio X tech repeats
                                [ibiorep,itechrep] = ind2sub(size(AnlOut.(iGene{1}).(iStrn{1}).(iCond{1}).(iSort{1}).FragLenHist.Reps),irep);
                                labelCollect = cat(2,labelCollect,[Exp '_' rawHist.Name num2str(ibiorep) '_' num2str(itechrep)]);
                                if ~isempty(FragLenHist.Reps{ibiorep,itechrep}) % because some reps might be empty
                                    plot([10:10:1000],FragLenHist.Reps{ibiorep,itechrep},'LineWidth',1.5,'Color',FragIn.Colors{icolor},'LineStyle',FragIn.CondStyles{icond});
                                    hold on;
                                end
                            end
                        end
                    end
                catch
                end
            end
        end
        
        title([Exp ' _ ' iGene{1}],'Interpreter','none');
        if FragIn.plotlegend==1
            legend(labelCollect,'FontSize',10,'Location','eastoutside','Interpreter','none');
        end
        
        if FragIn.plotlenVSmid == 0
            ylim([0 0.07]);
            xlim([0 1000]);
        end
        
        %% Drawing genomic features
        if FragIn.plotlenVSmid==1
            rectHeight=10;
            rectY=10;
            CDS(1) = min(AnlOut.CoordsTable{iGene{1},'CDS_5'},AnlOut.CoordsTable{iGene{1},'CDS_3'});
            CDS(2) = max(AnlOut.CoordsTable{iGene{1},'CDS_5'},AnlOut.CoordsTable{iGene{1},'CDS_3'});
            rectangle('Position',[max(CDS(1),StartStopPlot(1)),rectY,min(CDS(2),StartStopPlot(2)) - max(CDS(1),StartStopPlot(1)),rectHeight],'FaceColor',[0 0 0]);
%             annotation('textbox',[0.15 0.75 0.2 0.05],'String',strcat(AnlOut.CoordsTable{iGene{1},'CDS'},' CDS'),'Color',[0 0 0],'EdgeColor','none','FontSize',14);
            
            if ~isnan(AnlOut.CoordsTable{iGene{1},'U1_CDS_5'}) && ~isempty(AnlOut.CoordsTable{iGene{1},'U1_CDS_5'})
                U1_CDS(1) = min(AnlOut.CoordsTable{iGene{1},'U1_CDS_5'},AnlOut.CoordsTable{iGene{1},'U1_CDS_3'});
                U1_CDS(2) = max(AnlOut.CoordsTable{iGene{1},'U1_CDS_5'},AnlOut.CoordsTable{iGene{1},'U1_CDS_3'});
                try
                    rectangle('Position',[max(U1_CDS(1),StartStopPlot(1)),rectY+rectHeight/2,min(U1_CDS(2),StartStopPlot(2)) - max(U1_CDS(1),StartStopPlot(1)),rectHeight],'FaceColor',[0 0 1]);
%                     annotation('textbox',[0.15 0.8 0.2 0.05],'String',strcat(AnlOut.CoordsTable{iGene{1},'U1_CDS'},' CDS'),'Color',[0 0 1],'EdgeColor','none','FontSize',14);
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
                for iGal4 = 1:4
                    gal4bs(1) = min(Sites.Gal4bs.(iGene{1}){iGal4}(1),Sites.Gal4bs.(iGene{1}){iGal4}(2));
                    gal4bs(2) = max(Sites.Gal4bs.(iGene{1}){iGal4}(1),Sites.Gal4bs.(iGene{1}){iGal4}(2));
                    rectangle('Position',[gal4bs(1),rectY,gal4bs(2)-gal4bs(1)+1,rectHeight],'FaceColor',[0 0 1]);
                end
                for iMig1 = 1:2
                    mig1bs(1) = min(Sites.Mig1bs.(iGene{1}){iMig1}(1),Sites.Mig1bs.(iGene{1}){iMig1}(2));
                    mig1bs(2) = max(Sites.Mig1bs.(iGene{1}){iMig1}(1),Sites.Mig1bs.(iGene{1}){iMig1}(2));
                    rectangle('Position',[mig1bs(1),rectY,mig1bs(2)-mig1bs(1)+1,rectHeight],'FaceColor',[1 0 0]);
                end
                
                tata(1) = min(Sites.TATA.(iGene{1})(1),Sites.TATA.(iGene{1})(2));
                tata(2) = max(Sites.TATA.(iGene{1})(1),Sites.TATA.(iGene{1})(2));
                rectangle('Position',[tata(1),rectY,tata(2)-tata(1)+1,rectHeight],'FaceColor',[0 1 0]);
                
                tss(1) = min(Sites.TSS.(iGene{1})(1),Sites.TSS.(iGene{1})(2));
                tss(2) = max(Sites.TSS.(iGene{1})(1),Sites.TSS.(iGene{1})(2));
                rectangle('Position',[tss(1),rectY,tss(2)-tss(1)+1,rectHeight],'FaceColor',[1 0 1]);
            catch
            end
            
% configuring xlim            
            if U1_CDS(2)<CDS(1)
                xlim([max(U1_CDS(2)-200,StartStopPlot(1)), min(CDS(1)+200,StartStopPlot(2))]);
            elseif CDS(2)<U1_CDS(1)
                xlim([max(CDS(2)-200,StartStopPlot(1)), min(U1_CDS(1)+200,StartStopPlot(2))]);
                if strcmp(iGene{1},'HO_GAL1')
                    rectangle('Position',[CDS(2)+158,rectY-rectHeight/2,8,rectHeight*3],'FaceColor',[0.5 0 0]);
                    xlim([max(CDS(2)-200,StartStopPlot(1)), StartStopPlot(2)-500]); %StartStopPlot(2)-125
                    set(gca,'XDir','reverse');
                end
            end
            
            if FragIn.savefigs==1
                saveas(hfLenHist.(iGene{1}).(iStrn{1}),[FragIn.Dir '/LENvsMID_' iGene{1} '_' iStrn{1} '_' iCond{1} '_' iSort{1} '.png']);
                saveas(hfLenHist.(iGene{1}).(iStrn{1}),[FragIn.Dir '/LENvsMID_' iGene{1} '_' iStrn{1} '_' iCond{1} '_' iSort{1} '.fig']);
            end
        end
            
        if FragIn.plotlenVSmid==0 && FragIn.savefigs==1 && FragIn.stackStrns==0
            saveas(hfLenHist.(iGene{1}).(iStrn{1}),[FragIn.Dir '/LenHist_' iGene{1} '_' iStrn{1} '.png']);
            saveas(hfLenHist.(iGene{1}).(iStrn{1}),[FragIn.Dir '/LenHist_' iGene{1} '_' iStrn{1} '.fig']);
        end
    end
    if FragIn.plotlenVSmid==0 && FragIn.savefigs==1 && FragIn.stackStrns==1
        saveas(hfLenHist.(iGene{1}),[FragIn.Dir '/LenHist_' iGene{1} '.png']);
        saveas(hfLenHist.(iGene{1}),[FragIn.Dir '/LenHist_' iGene{1} '.fig']);
    end
end
if FragIn.savefigs==1
    CurrentCodeFile = [mfilename '.m'];
    [status,msg] = copyfile(CurrentCodeFile,FragIn.Dir);
    if status==0
        error('Error. hist-plotting code file was not copied')
    end
    save([FragIn.Dir '/FragIn_' datestr(now,'yymmdd_HHMMSS') '.mat'],'FragIn');
end