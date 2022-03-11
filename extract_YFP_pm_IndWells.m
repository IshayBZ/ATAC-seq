function logYFP = extract_YFP_pm_IndWells(YFPIn)
%% Extract ON fraction - ON level metrics per Stratedigm run (20200221)
% Extract individual wells as they're arranged in the plate (not in the flipped DG format), plot, and save in struct

if YFPIn.saveoutput==1
    outdir = [YFPIn.datadir 'out_' datestr(now,'yymmdd_HHMM') '/'];
    mkdir(outdir);
end

if YFPIn.filDeb==1
    load(YFPIn.DebGates);
else
    gates=[];
end

if YFPIn.plot_means==1
    hfMeans = figure('WindowStyle','docked');
    tile = tiledlayout('flow');
    tile.TileSpacing = 'compact';
    tile.Padding = 'compact';
    axMeans = nexttile(tile);
end
%% Gauss plot params
cm = parula; m0 = 0; mmax = 1;

RsqTh = 0.99; %Threshold for fitting two Gaussians and for black title
RsqTh2 = 0.9; %Second threshold (only for ploting)

par.muL = 0.2; par.muS = 2; par.muU = 4;
par.sigL = 0.05; par.sigS = 0.13; par.sigU = 0.5;

par.alphaL = 0; par.alphaS = 0.5; par.alphaU = 1;
par.mu1L = 0.2; par.mu1S = 0.5; par.mu1U = 1.2;
par.sig1L = 0.05; par.sig1S = 0.13; par.sig1U = 0.5;
par.mu2L = 1.2; par.mu2S = 2.5; par.mu2U = 4;
par.sig2L = 0.05; par.sig2S = 0.2; par.sig2U = 0.5;
if isfield(YFPIn,'muTH')
    par.mu1U = YFPIn.muTH;
    par.mu2L = YFPIn.muTH;
end
    
%%
yfpedges = linspace(0,4,101);
yfprange = yfpedges(1:end-1);

StrntoGene = readtable([YFPIn.datadir YFPIn.PlateMap],'ReadRowNames',true,'Sheet','StrntoGene');
for platename = YFPIn.Plates
    if YFPIn.plot_hist==1
        [hfHist ha] = gridplot_docked(8,12,60,60,'gap',2);
    end
    PlateMap = readtable([YFPIn.datadir YFPIn.PlateMap],'ReadRowNames',true,'Sheet', platename{1});
    
    %% Plotting and extracting
    % UPDATED 200226 IBZ - changed the order of extraction and plotting to fit the actual plate layout (and not flip vertically as in DG) but meaning of r,c is unchanged
    for r = 1:8 %8:-1:rmin
        for c = 1:12
            if YFPIn.plot_hist==1
                if isfield(YFPIn,'DGformat') % i.e. YFPIn.DGformat = 1
                    axes(ha(sub2ind([12 8],c,9-r))); % to plot DGs as the flipped lab's format
                else
                    axes(ha(sub2ind([12 8],c,r)));
                end
            end
            [data,meta] = fcsparsewell(YFPIn.datadir, platename{1}, [r c]);
            
            %check for empty/noisy (low OD) wells
            if isempty(data) || fcsnumel(data)==0
                continue;
            end
            
            %filter stage 1 (discard off-scale)
            datafilt = fcsselect(data, data.fsc<9999 & data.ssc < 9999);
            %filter stage 2: manual debris gating
            if YFPIn.filDeb==1
                subpops = applygate(datafilt,gates);
                datafiltgated = subpops(2);
            else
                datafiltgated = datafilt; %ignore filter stage 2
            end
            
            logyfpdat = log10((10^3.5)*(datafiltgated.yfp./datafiltgated.ssc));
            pm = mean(logyfpdat);
            
            %Creating a PDF version of yfphist
            [yfppdf,~] = histcounts(logyfpdat,yfpedges,'Normalization','pdf');
            if YFPIn.plot_hist==1
                plot(yfprange,yfppdf,'k-','linewidth',2); hold on;
            end
            
            %%% samples missing from platemap were plotted but will not be saved and won't display title and means
            if strcmp(PlateMap{r,c}{1},'-')
                continue
            end
            
            smplName = regexp(PlateMap{r,c},'_','split'); %Names in this platemap shouldn't contain Tagmentation or sorting information
            iStrn = [smplName{1}{1} '_' smplName{1}{2}];
            if strncmp(iStrn,'S_GAL',5) % for GAL reporters YFP data that I'm corelating with WT ATAC-seq data
                iCorelStrn = 'S_W';
            else
                iCorelStrn = iStrn;
            end
            iCond = smplName{1}{3};
            iBioRep = str2num(smplName{1}{4});
            %                 if ~isempty(str2num(smplName{1}{end-1}(end))) %meaning it's iBioRep and it's followed by iTechRep
            %                     iBioRep = str2num(smplName{1}{end-1}(end));
            %                 else %meaning there's no iTechRep
            %                     iBioRep = str2num(smplName{1}{end}(end));
            %                 end
            iGene = StrntoGene.Reporter(iStrn);
            logYFP.(iGene{1}).(iCorelStrn).(iCond).PM{iBioRep} = pm; % Population mean
            logYFP.(iGene{1}).(iCorelStrn).(iCond).DATA{iBioRep} = logyfpdat;
            
            %% Fitting one and two Gaussians
            % checking if sample should be bimodal
            bimConds = {'D0Gp1','Dm2Gp1','Dm1G0','Dm2Gm1','Dm3Gm2','Dm2Gm2','Dm2G0','Dm4Gm3'};
            if (strncmp(iCorelStrn,'S_W',3) || strcmp(iCorelStrn,'S_M')) && ~isequal(strcmp(iCond,bimConds),zeros(1,length(bimConds)))
                % i.e. if iStrn = S_W or S_WL93 etc && iCond = one of bimConds
                Gauss1 = @(mu,sig,x) exp(-(((x-mu).^2)/(2*sig^2)))/sqrt(2*pi*sig^2);
                Gauss2 = @(alpha,mu1,sig1,mu2,sig2,x) ...
                    alpha*exp(-(((x-mu1).^2)/(2*sig1^2)))/sqrt(2*pi*sig1^2) +...
                    (1-alpha)*exp(-(((x-mu2).^2)/(2*sig2^2)))/sqrt(2*pi*sig2^2);
                
                [fitobject,gof,~] = fit(yfprange',yfppdf',Gauss1,'Lower',[par.muL par.sigL],...
                    'StartPoint',[par.muS par.sigS],'Upper',[par.muU par.sigU]);
                
                if gof.adjrsquare >= RsqTh %Sticking with one Gaussian fit
                    if fitobject.mu < 1
                        of = 0;
                        fm = fitobject.mu;
                        fs = fitobject.sig;
                        om = NaN;
                        os = NaN;
                    else
                        of = 1;
                        fm = NaN;
                        fs = NaN;
                        om = fitobject.mu;
                        os = fitobject.sig;
                    end
                    if YFPIn.print_title==1 && YFPIn.plot_hist==1
                        title(['mu = ',num2str(fitobject.mu,2)],'fontsize',10,'Color','k','VerticalAlignment','top');
                    end
                elseif gof.adjrsquare < RsqTh % then fit two Gaussians
                    [fitobject,gof,~] = fit(yfprange',yfppdf',Gauss2,...
                        'Lower',[par.alphaL par.mu1L par.sig1L par.mu2L par.sig2L],...
                        'StartPoint',[par.alphaS par.mu1S par.sig1S par.mu2S par.sig2S],...
                        'Upper',[par.alphaU par.mu1U par.sig1U par.mu2U par.sig2U]);
                    of = 1-fitobject.alpha;
                    fm = fitobject.mu1;
                    fs = fitobject.sig1;
                    om = fitobject.mu2;
                    os = fitobject.sig2;
                    if YFPIn.print_title==1 && YFPIn.plot_hist==1
                        if gof.adjrsquare >= RsqTh %(of the two Gaussians fit)
                            title([num2str(fm,2),' ',num2str(om,2)],'fontsize',10,'Color','k','VerticalAlignment','top');
                        elseif gof.adjrsquare < RsqTh && gof.adjrsquare >= RsqTh2
                            title([num2str(fm,2),' ',num2str(om,2)],'fontsize',10,'Color','m','VerticalAlignment','top');
                        elseif gof.adjrsquare < RsqTh2
                            title([num2str(fm,2),' ',num2str(om,2),'   R^2=',num2str(gof.adjrsquare,2)],...
                                'fontsize',10,'Color','r','VerticalAlignment','top');
                        end
                        if fm > om
                            title('\mu_1 > \mu_2','fontsize',12);
                        end
                    end
                end
                if YFPIn.print_title==1 && YFPIn.plot_hist==1
                    plot(fitobject,yfprange',yfppdf');%,'b-','linewidth',2);
                    %                     hold on
                    plot([om om],[0 5],'Color',[0.5 0.5 0.5]);
                    axHist = gca;
                    axHist.Color = cm(min(length(cm),floor((of-m0)/(mmax-m0)*length(cm))+1),:);
                end
                logYFP.(iGene{1}).(iCorelStrn).(iCond).OF{iBioRep} = of; % ON fraction
                logYFP.(iGene{1}).(iCorelStrn).(iCond).FM{iBioRep} = fm; % OFF mean
                logYFP.(iGene{1}).(iCorelStrn).(iCond).FS{iBioRep} = fs; % OFF sigma
                logYFP.(iGene{1}).(iCorelStrn).(iCond).OM{iBioRep} = om; % ON mean
                logYFP.(iGene{1}).(iCorelStrn).(iCond).OS{iBioRep} = os; % ON sigma
            else
                if YFPIn.print_title==1 && YFPIn.plot_hist==1
                    title(num2str(pm,2),'fontsize',8,'Color','k','VerticalAlignment','top');
                end
            end
            
            if YFPIn.plot_hist==1
                plot([pm pm],[0 5],'Color','r');
                legend('off');
                xlabel('')
                ylabel('')
                xlim([0 4]);
                ylim([0 5]);
                axHist = gca;
                axHist.XTickLabel = [];
                axHist.YTickLabel = [];
            end
            
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
                if isfield(logYFP.(iGene{1}).(iCorelStrn).(iCond),'OF')
                    plot(Xjitt,fm,'Color',YFPIn.Colors.(iCorelStrn),'Marker',YFPIn.Markers{iBioRep},'MarkerSize',5,'LineWidth',2); hold on;
                    plot(Xjitt,om,'Color',YFPIn.Colors.(iCorelStrn),'Marker',YFPIn.Markers{iBioRep},'MarkerSize',10,'LineWidth',2);
                else
                    plot(Xjitt,pm,'Color',YFPIn.Colors.(iCorelStrn),'Marker',YFPIn.Markers{iBioRep},'MarkerSize',10,'LineWidth',2);
                end
                hold on;
            end
        end
    end
    
    if YFPIn.plot_hist==1 && YFPIn.saveoutput==1
        saveas(hfHist,[outdir platename{1} '_Hist.fig']);
        saveas(hfHist,[outdir platename{1} '_Hist.png']);
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
    save([outdir '/extracted metrics'],'logYFP','gates');
    % save code snapshot in the same folder
    CurrentCodeFile = [mfilename '.m'];
    [status,msg] = copyfile(CurrentCodeFile,outdir);
    if status==0
        error('Error. code file was not copied')
    end
end