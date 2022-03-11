function logYFP = extract_YFP_from_old_metrics(YFPIn)
%% Extract ON fraction - ON level metrics per Stratedigm run (20200221)
% Extract individual wells as they're arranged in the plate (not in the flipped DG format), plot, and save in struct

if YFPIn.saveoutput==1
    outdir = [YFPIn.datadir 'out_' datestr(now,'yymmdd_HHMM') '/'];
    mkdir(outdir);
end

if YFPIn.plot_means==1
    hfMeans = figure('WindowStyle','docked');
    tile = tiledlayout('flow');
    tile.TileSpacing = 'compact';
    tile.Padding = 'compact';
    axMeans = nexttile(tile);
end
%%
StrntoGene = readtable([YFPIn.datadir YFPIn.PlateMap],'ReadRowNames',true,'Sheet','StrntoGene');
istrn = 0;
for iStrn = StrntoGene.Properties.RowNames'
    istrn = istrn + 1;
    iGene = StrntoGene.Reporter(iStrn{1});
    Plates = regexp(StrntoGene.Plates(iStrn{1}),',','split');
    OldExtMetrics = StrntoGene.OldExtMetrics(iStrn{1});
    if strncmp(iStrn{1},'S_GAL',5) % for GAL reporters YFP data that I'm corelating with WT ATAC-seq data
        iCorelStrn = {'S_W'};
    else
        iCorelStrn = iStrn;
    end
    load([YFPIn.OldExtPath OldExtMetrics{1}]);
    for platename = Plates{1}
        PlateMap = readtable([YFPIn.datadir YFPIn.PlateMap],'ReadRowNames',true,'Sheet', platename{1});
        
        for r = 8:-1:1 %going from last to first (cause format of DGs was flipped in analysis)
            for c = 1:12
                if strcmp(PlateMap{r,c}{1},'-') % platemaps reflect the actual (experimental) layout
                    continue
                end
                
                smplName = regexp(PlateMap{r,c},'_','split'); %Names in this platemap shouldn't contain Tagmentation or sorting information
                iCond = smplName{1}{3};
                iBioRep = str2num(smplName{1}{4});
                
                % using 9-r since rows were flipped in these arrays to resemble the DG figure layout
                % but my plate layout is not flipped and resembles the experimental layout
                logYFP.(iGene{1}).(iCorelStrn{1}).(iCond).PM{iBioRep} = collectPM(istrn,iBioRep,9-r,c); % Population mean
                logYFP.(iGene{1}).(iCorelStrn{1}).(iCond).DATA{iBioRep} = alldata(istrn,iBioRep,9-r,c); % all YFP data
                logYFP.(iGene{1}).(iCorelStrn{1}).(iCond).OF{iBioRep} = collectOF(istrn,iBioRep,9-r,c); % ON fraction
                logYFP.(iGene{1}).(iCorelStrn{1}).(iCond).FM{iBioRep} = collectFM(istrn,iBioRep,9-r,c); % OFF mean
                logYFP.(iGene{1}).(iCorelStrn{1}).(iCond).FS{iBioRep} = collectFS(istrn,iBioRep,9-r,c); % OFF sigma
                logYFP.(iGene{1}).(iCorelStrn{1}).(iCond).OM{iBioRep} = collectOM(istrn,iBioRep,9-r,c); % ON mean
                logYFP.(iGene{1}).(iCorelStrn{1}).(iCond).OS{iBioRep} = collectOS(istrn,iBioRep,9-r,c); % ON sigma
                
                if YFPIn.plot_means==1
                    if YFPIn.plotgrad==1 % gluc/gal gradient
                        % [~,~,Gluc,Gal,~] = nameConvert_ATAC_to_DG(iCond);
                        [~,Gluc,~] = nameConvert_ATAC_to_DG_210602(iCond);
                        %                         X = Gluc;
                        %                         XjitSize = 0.01 + X/20;
                        X = max(-6,log2(Gluc));
                        XjitSize = 0.5 + X/20;
                    elseif YFPIn.plotgrad==2 % estradiol gradient
                        [~,Est,~,~] = nameConvert_ATAC_to_DG_211108(iCond);
                        X = Est;
                        XjitSize = 6 + X/20;
                    elseif YFPIn.plotgrad==3 % dox gradient
                        [Dox,~,~,~] = nameConvert_ATAC_to_DG_211108(iCond);
                        X = max(0,log2(Dox));
                        XjitSize = 0.5 + X/20;
                    end
                    
                    Xjitt = X + XjitSize * rand(1,length(X)) - XjitSize/2;
                    
                    axes(axMeans);
                    if isfield(logYFP.(iGene{1}).(iCorelStrn{1}).(iCond),'OF') % NEED TO CORRECT - OF is always a field but might be empty
                        plot(Xjitt,fm,'Color',YFPIn.Colors.(iCorelStrn{1}),'Marker',YFPIn.Markers{iBioRep},'MarkerSize',5,'LineWidth',2); hold on;
                        plot(Xjitt,om,'Color',YFPIn.Colors.(iCorelStrn{1}),'Marker',YFPIn.Markers{iBioRep},'MarkerSize',10,'LineWidth',2);
                    else
                        plot(Xjitt,pm,'Color',YFPIn.Colors.(iCorelStrn{1}),'Marker',YFPIn.Markers{iBioRep},'MarkerSize',10,'LineWidth',2);
                    end
                    hold on;
                end
            end
        end
    end
end
if YFPIn.plot_means==1
    if YFPIn.plotgrad==1
        xlabel(axMeans,'max(-6,log2(Gluc))');
    elseif YFPIn.plotgrad==2
        xlabel(axMeans,'Est');
    elseif YFPIn.plotgrad==3
        xlabel(axMeans,'max(0,log2(Dox))');
    end
    ylabel(axMeans,'log10(YFP)');
    
    if YFPIn.saveoutput==1
        saveas(hfMeans,[outdir 'Means.fig']);
        saveas(hfMeans,[outdir 'Means.png']);
    end
end

if YFPIn.saveoutput==1
    save([outdir '/extracted metrics'],'logYFP');
    % save code snapshot in the same folder
    CurrentCodeFile = [mfilename '.m'];
    [status,msg] = copyfile(CurrentCodeFile,outdir);
    if status==0
        error('Error. code file was not copied')
    end
end