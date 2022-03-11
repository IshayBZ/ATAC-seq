function plot_ATAC_YFP_corel(Exp,YFPExp,CorelIn,logYFP,AnlOut,Regs)

for iGene = CorelIn.Genes % For GAL10 use GAL1 ATACseq data and GAL10 expression data
    hfCorel.(iGene{1}) = figure('WindowStyle','docked');
    tile.(iGene{1}) = tiledlayout(2,4); %tiledlayout('flow');
    tile.(iGene{1}).TileSpacing = 'compact';
    tile.(iGene{1}).Padding = 'compact';
    labelCollect={}; % collect labels
    
    if strcmp(iGene{1},'GAL10')
        iATACgene = 'GAL1'; %if looking at GAL10 expression use GAL1 for accessibility
    else
        iATACgene = iGene{1};
    end
    
    meanlogTable.(iGene{1}) = table; % this table is made once but overwritten (with identical entries) for every iReg (since I'm plotting a panel for each iReg)
    for iReg = fieldnames(Regs.(iGene{1}))'
        reg = Regs.(iGene{1}).(iReg{1}) - AnlOut.(iATACgene).StartStop(1);
        ax = nexttile(tile.(iGene{1}));
        
        if strcmp(iGene{1},'HO_GAL1')
            iYFPgene = 'GAL1'; %if looking at HO_GAL1 for accessibility, use GAL1 for expression
            rawHistON = AnlOut.(iATACgene).rawHistONSNP.Avg(reg(1):reg(2));
            rawHistOFF = AnlOut.(iATACgene).rawHistOFFSNP.Avg(reg(1):reg(2));
        else
            iYFPgene = iGene{1};
            rawHistON = AnlOut.(iATACgene).rawHistON.Avg(reg(1):reg(2));
            rawHistOFF = AnlOut.(iATACgene).rawHistOFF.Avg(reg(1):reg(2));
        end
        
        sumhistTable.(iGene{1}).(iReg{1}) = table;
        sumhistscaledTable.(iGene{1}).(iReg{1}) = table;
        RsqCollect = [];
        
        if CorelIn.ConcatStrns==1 % concatenates strns into one (e.g. from diff Xlink experiments such as S_W and S_WL93)
            istrain = ['S_' CorelIn.Strns{1}]; % takes the first entry in CorelIn.Strns
            meanlog.(iGene{1}).(istrain).nonzero = []; % this is calculated all over again for each iReg
            sumhist.(iGene{1}).(iReg{1}).(istrain).nonzero = [];
            sumhistscaled.(iGene{1}).(iReg{1}).(istrain).nonzero = [];
        end
        for iStrn = CorelIn.Strns %fieldnames(logYFP.(iGene{1}))'
            iStrn = {['S_' iStrn{1}]};
            if CorelIn.ConcatStrns==0 % if not concatenating strains, save each strain in its own struct field
                istrain = iStrn{1};
                meanlog.(iGene{1}).(istrain).nonzero = []; % this is calculated all over again for each iReg
                sumhist.(iGene{1}).(iReg{1}).(istrain).nonzero = [];
                sumhistscaled.(iGene{1}).(iReg{1}).(istrain).nonzero = [];
            end
            if strcmp(iStrn{1},'S_WL93') % since YFP data is always saved as S_W
                iYFPStrn = {'S_W'};
            else
                iYFPStrn = iStrn;
            end
            for iCond = CorelIn.Conds %fieldnames(logYFP.(iGene{1}).(iStrn{1}))'
                if strncmp(iStrn{1},'S_W',3) && strcmp(iCond{1},'Dp1') && strncmp(YFPExp,'L13',3)
                    iYFPCond = {'D0'}; %Dp1 doesn't exist in DGs, so use D0 which is a good enough estimate
                else
                    iYFPCond = iCond;
                end
                for iSort = CorelIn.Sorts
                    MEANLOG = zeros(4,8);
                    SUMHIST = zeros(4,8);
                    for iBioRep = 1:4
                        if isfield(CorelIn,'TechReps')
                            TechReps = CorelIn.TechReps;
                        else
                            TechReps = 1:8;
                        end
                        for iTechRep = TechReps % note that all have the same YFP value
                            if iTechRep==1
                                MarkFaceCol = [1 1 1];
                            else
                                MarkFaceCol = CorelIn.Colors.(iStrn{1});
                            end
                            try % empty bioreps will be skipped, as well as samples with no YFP data
                                if isfield(logYFP.(iYFPgene).(iYFPStrn{1}).(iYFPCond{1}),'FM')
                                    if strcmp(iSort,'F')
                                        MEANLOG(iBioRep,iTechRep) = logYFP.(iYFPgene).(iYFPStrn{1}).(iYFPCond{1}).FM{iBioRep};
                                    elseif strcmp(iSort,'O')
                                        MEANLOG(iBioRep,iTechRep) = logYFP.(iYFPgene).(iYFPStrn{1}).(iYFPCond{1}).OM{iBioRep};
                                    end
                                    if isnan(MEANLOG(iBioRep,iTechRep)) %if a supposedly bimodal sample is actually unimodal then sort F will be NaN --> ignore this sort
                                        MEANLOG(iBioRep,iTechRep) = 0;
                                    end
                                else
                                    MEANLOG(iBioRep,iTechRep) = logYFP.(iYFPgene).(iYFPStrn{1}).(iYFPCond{1}).PM{iBioRep};
                                end
                            catch end
                            try % empty reps will be skipped, as well as samples with no ATAC data
                                rawHist = AnlOut.(iATACgene).(iStrn{1}).(iCond{1}).(iSort{1}).rawHist.Reps{iBioRep,iTechRep}(reg(1):reg(2));
                                if CorelIn.ScaleATAC==0
                                    SUMHIST(iBioRep,iTechRep) = sum(rawHist)/(reg(2) - reg(1) + 1);
                                    fitslopestartpoint = 10^4;
                                elseif CorelIn.ScaleATAC==1
                                    SUMHIST(iBioRep,iTechRep) = sum(rawHist) ./ sum(rawHistON);
                                    fitslopestartpoint = 5;
                                end
                            catch end
                            %plotting direct correlations of corresponding reps
                            if MEANLOG(iBioRep,iTechRep)*SUMHIST(iBioRep,iTechRep)~=0 && CorelIn.ErrBar==0 && strcmp(Exp,YFPExp)
                                plot(SUMHIST(iBioRep,iTechRep),10^MEANLOG(iBioRep,iTechRep),'Color',CorelIn.Colors.(iStrn{1}),...
                                    'Marker',CorelIn.Markers{iBioRep},'MarkerFaceColor',MarkFaceCol,'MarkerSize',10,'LineWidth',2);
                                hold on;
                            end
                            %gather nonzero ATAC data in a table
                            if SUMHIST(iBioRep,iTechRep)~=0
                                smplname = [iCond{1} '_' iSort{1} num2str(iBioRep) '_' num2str(iTechRep)];
                                sumhistTable.(iGene{1}).(iReg{1}){smplname,istrain} = SUMHIST(iBioRep,iTechRep);
                            end
                        end
                        %gather nonzero YFP data in a table
                        if MEANLOG(iBioRep,iTechRep)~=0 %this is the last iTechRep in the loop
                            smplname = [iCond{1} '_' iSort{1} num2str(iBioRep)];
                            meanlogTable.(iGene{1}){smplname,istrain} = MEANLOG(iBioRep,iTechRep); % this value is overwritten (with identical entries) again and again for each iReg
                        end
                    end
                    % gather nonzero points in variables for fitting direct correlations
                    meanlog.(iGene{1}).(istrain).nonzero = [meanlog.(iGene{1}).(istrain).nonzero; MEANLOG(find(MEANLOG.*SUMHIST~=0))];
                    sumhist.(iGene{1}).(iReg{1}).(istrain).nonzero = [sumhist.(iGene{1}).(iReg{1}).(istrain).nonzero; SUMHIST(find(MEANLOG.*SUMHIST~=0))];
                    % gather arrays in variables
                    meanlog.(iGene{1}).(istrain).(iCond{1}).(iSort{1}) = MEANLOG;
                    sumhist.(iGene{1}).(iReg{1}).(istrain).(iCond{1}).(iSort{1}) = SUMHIST;
                    % plot averages + error bars over all reps
                    if CorelIn.ErrBar==1
                        SUMHISTavg = mean(nonzeros(SUMHIST));
                        SUMHISTerr = std(nonzeros(SUMHIST)) / sqrt(length(nonzeros(SUMHIST)));
                        MEANLOGavg = mean(nonzeros(MEANLOG(:,iTechRep))); %iTechRep is the last rep in the loop
                        MEANLOGerr = std(nonzeros(MEANLOG(:,iTechRep))) / sqrt(length(nonzeros(MEANLOG(:,iTechRep))));
                        % errorbar(x,y,yneg,ypos,xneg,xpos)
%                         errorbar(SUMHISTavg,MEANLOGavg,MEANLOGerr,MEANLOGerr,SUMHISTerr,SUMHISTerr,'o');
                        yneg = 10^MEANLOGavg * (1 - 1/(10^MEANLOGerr));
                        ypos = 10^MEANLOGavg * (10^MEANLOGerr - 1);
                        errorbar(SUMHISTavg,10^MEANLOGavg,yneg,ypos,SUMHISTerr,SUMHISTerr,'Color',CorelIn.Colors.(iStrn{1}),'Marker','o','LineWidth',1);
                        hold on;
                    end
                end
            end
            % fit in lin or log scale - CORRECT TO ALLOW FITTING OF AVERAGES FOR INDIRECT CORRELATIONS
            fitfun = 'A*x+B';
            if CorelIn.ConcatStrns==0
                SUMHISTFIT = linspace(min(sumhist.(iGene{1}).(iReg{1}).(istrain).nonzero),max(sumhist.(iGene{1}).(iReg{1}).(istrain).nonzero),100);
                if strcmp(CorelIn.YFPscale,'log')
                    [fobj,gof] = fit(sumhist.(iGene{1}).(iReg{1}).(istrain).nonzero,meanlog.(iGene{1}).(istrain).nonzero,fitfun,'StartPoint',[fitslopestartpoint 0.5]);
                    MEANLOGFIT = fobj(SUMHISTFIT);
                    plot(SUMHISTFIT,10.^MEANLOGFIT,'-','linewidth',1,'Color',CorelIn.Colors.(iStrn{1}));
                    RsqCollect = [RsqCollect '   R^2_' (istrain(3:end)) ' = ' num2str(gof.adjrsquare,2)];
                elseif strcmp(CorelIn.YFPscale,'linear')
                    [fobj,gof] = fit(sumhist.(iGene{1}).(iReg{1}).(istrain).nonzero,10.^meanlog.(iGene{1}).(istrain).nonzero,fitfun,'StartPoint',[10*fitslopestartpoint 0]);
                    MEANFIT = fobj(SUMHISTFIT);
                    plot(SUMHISTFIT,MEANFIT,'-','linewidth',1,'Color',CorelIn.Colors.(iStrn{1}));
                    RsqCollect = [RsqCollect '   R^2_' (istrain(3:end)) ' = ' num2str(gof.adjrsquare,2)];
                end
            end
        end
        if CorelIn.ConcatStrns==1
            SUMHISTFIT = linspace(min(sumhist.(iGene{1}).(iReg{1}).(istrain).nonzero),max(sumhist.(iGene{1}).(iReg{1}).(istrain).nonzero),100);
            if strcmp(CorelIn.YFPscale,'log')
                [fobj,gof] = fit(sumhist.(iGene{1}).(iReg{1}).(istrain).nonzero,meanlog.(iGene{1}).(istrain).nonzero,fitfun,'StartPoint',[fitslopestartpoint 0.5]);
                MEANLOGFIT = fobj(SUMHISTFIT);
                plot(SUMHISTFIT,10.^MEANLOGFIT,'-','linewidth',1,'Color',CorelIn.Colors.(istrain));
                RsqCollect = [RsqCollect '   R^2_' (istrain(3:end)) ' = ' num2str(gof.adjrsquare,2)];
            elseif strcmp(CorelIn.YFPscale,'linear')
                [fobj,gof] = fit(sumhist.(iGene{1}).(iReg{1}).(istrain).nonzero,10.^meanlog.(iGene{1}).(istrain).nonzero,fitfun,'StartPoint',[10*fitslopestartpoint 0]);
                MEANFIT = fobj(SUMHISTFIT);
                plot(SUMHISTFIT,MEANFIT,'-','linewidth',1,'Color',CorelIn.Colors.(istrain));
                RsqCollect = [RsqCollect '   R^2_' (istrain(3:end)) ' = ' num2str(gof.adjrsquare,2)];
            end
        end
        title([iGene{1} ' - ' iReg{1} RsqCollect],'Interpreter','none','FontSize',8);
        %         legend(labelCollect,'FontSize',10,'Location','eastoutside','Interpreter','none');
        %         legend('off'); % I made the legened to have data points names in the plot browser but I don't want it to appear all the time
        
        xlabel('Accessibility'); %,'FontSize',20);
        if CorelIn.ScaleATAC==1
            xlim([0 1.2]);
        else
%             xlim([4E-4 18E-4]);
        end
        ylabel('YFP'); %,'FontSize',20);
        ylim([3 10000]);
        
        ax.YScale = CorelIn.YFPscale;
        
        %save tables to excel
        if CorelIn.savefigs==1
            writetable(sumhistTable.(iGene{1}).(iReg{1}),[CorelIn.Dir '/Table.xlsx'],'WriteRowNames',true,'Sheet', [(iGene{1}) '_Acc_' (iReg{1})]);
            %if this doesnt work make a local dir and print this note 'PATH TOO LONG - COREL FIGS SAVED IN CURRENT FOLDER'
        end
    end
    
    if CorelIn.savefigs==1
        writetable(meanlogTable.(iGene{1}),[CorelIn.Dir '/Table.xlsx'],'WriteRowNames',true,'Sheet', [iGene{1} '_log10(YFP)']);
        saveas(hfCorel.(iGene{1}),[CorelIn.Dir '/Corel_' iGene{1} '.png']);
        saveas(hfCorel.(iGene{1}),[CorelIn.Dir '/Corel_' iGene{1} '.fig']);
        %     close(gcf);
    end
end
if CorelIn.savefigs==1
    save([CorelIn.Dir '/CorelIn_' datestr(now,'yymmdd_HHMMSS') '.mat'],'CorelIn','sumhist','sumhistTable','meanlog','meanlogTable');
    CurrentCodeFile = [mfilename '.m'];
    [status,msg] = copyfile(CurrentCodeFile,CorelIn.Dir);
    if status==0
        error('Error. plotting code file was not copied')
    end
end