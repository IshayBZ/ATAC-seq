%% ATACseq analysis pipeline
clearvars
% PURPOSE
% 1. EXTRACT ATAC-SEQ INSERTIONS AS GIVEN BY FRAGMENT ENDS IN TXT FILES
% 2. ARRANGE IN STRUCTS OF gene X strain X condition X sort gate
% 3. PLOT INSERTION HISTOGRAMS & SUMS OVER SPECIFIED REGIONS
% 4. EXTRACT & PLOT FLOW CYTOMETRY DATA & PLOT CORRELATION WITH ATAC-SEQ SUMS

% TO DO
% 1. CHOOSE ATAC-SEQ AND FLOW CYTOMETRY EXPERIMENT NUMBERS
% 2. SPECIFY PARAMS IN RELEVANT SECTION
% 3. SPECIFY PATHS OF COMPLETED STEPS TO SKIP AND SAVE TIME

Exp = 'O16'; YFPExp = 'O16';

%% L93 params
if strcmp(Exp,'L93')
    datadir = '../..\O12 - L93 Reanalysis\data/';
    % NOTE: rerun data is not analyzed (relevant to gal80del and double-del)
    SmplNameON = 'S_W_Gp1_A'; % Used to calculate the scaled sum 'sumHistScaled'
    SmplNameOFF = 'S_W_Dp1_A'; % Used for both 'sumHistScaled' AND to plot histograms divided by the OFF histogram
    SmplNameONSNP = SmplNameON; %work around since there were no SNPs in this Exp
    SmplNameOFFSNP = SmplNameOFF;
    AnlIn.ReadThresh = 10000;
    AnlIn.BadSmpls = {};
    AnlIn.ut = 'uniq'; % 'uniq'/'tot' = remove/don't remove duplicates
    AnlIn.Norm = 'Loci'; % 'Loci'/'Chr' = Normalize by control loci / uncaptured chromosomes
    ExtDir = [datadir 'Ext_220204_1701_NF'];
%       AnlDir = [ExtDir '/Anl_220204_1702_uniq_NormChr_Th10000'];
    AnlDir = [ExtDir '/Anl_220204_1706_uniq_NormLoci_Th10000'];
    
    PlotHists = 0; PlotFrags = 0;
    HistIn.Genes = {'GAL1'}; HistIn.Strns = {'TETM'};
    HistIn.Conds = {'X0','X10','X20','X39','X78','X156','X313','X625','X1250','X2500','X5000','X10000'}; HistIn.Sorts = {'A'};
    HistIn.ColorBySort = 0; HistIn.CondStyles = {'-','-','-','-','-','-','-','-'};
    HistIn.ColorByCond = 1; HistIn.SortStyles = {'-','--','.-'};
    
    PlotSums = 0; SumIn.Genes = HistIn.Genes; SumIn.Strns = HistIn.Strns; SumIn.Conds = HistIn.Conds; SumIn.Sorts = HistIn.Sorts;
    PlotCorels = 1;
    % Tet-Mig1 correlation with SW25-26 YFP data added from O19
    yfpOutDir = '../../O12 - L93 Reanalysis\Stratedigm\Tet-MIG1 and SW25-26\out_220204_1647_FITC2';
        YFPIn.datadir = '../../O12 - L93 Reanalysis\Stratedigm\Tet-MIG1 and SW25-26\';
        YFPIn.PlateMap = 'platemap_220204.xlsx';
        YFPIn.Plates = {'plate1','plate'}; % plate1 = 2% FITC, plate2 = 4% FITC gain
        YFPIn.Genes = {'GAL1'};
        YFPIn.plot_hist=1; YFPIn.print_title=1; YFPIn.plot_means=1; YFPIn.plotgrad = 3; % choose gradient for PM plot 1=gluc/gal, 2=estradiol, 3=dox
    CorelIn.Genes = {'GAL1'}; CorelIn.Strns = {'TETM','W'}; CorelIn.ConcatStrns = 0; % concatenates strns into one (e.g. from diff Xlink experiments such as W and WL93)
    CorelIn.Conds = {'X0','X10','X20','X39','X78','X156','X313','X625','X1250','X2500','X5000','X10000','Dp1','Gp1'}; CorelIn.Sorts = {'A'};
    CorelIn.ScaleATAC = 1; CorelIn.ErrBar = 1; CorelIn.YFPscale = 'linear'; % set 'log' / 'linear' YScale
end
%% L108 params
if strcmp(Exp,'L108')
    datadir = '../..\O25 - L108 Reanalysis\data/align_220126/';
    SmplNameON = 'S_WL93_Gp1_O'; % Used to calculate the scaled sum 'sumHistScaled'
    SmplNameOFF = 'S_WL93_Dp1_F'; % Used for both 'sumHistScaled' AND to plot histograms divided by the OFF histogram
    SmplNameONSNP = SmplNameON; %work around since there were no SNPs in this Exp
    SmplNameOFFSNP = SmplNameOFF;
    AnlIn.ReadThresh = 5000;
    AnlIn.BadSmpls = {};
    AnlIn.ut = 'uniq'; % 'uniq'/'tot' = remove/don't remove duplicates
    AnlIn.Norm = 'Loci'; % 'Loci'/'Chr'/'Uncap' = Normalize by control loci / chromosomes 3,11,14 / all uncaptured reads
    ExtDir = [datadir 'Ext_220131_2240_NF']; %both normalizations are saved from now on
        AnlDir = [ExtDir '/Anl_220201_0941_uniq_NormLoci_Th5000'];
%         AnlDir = [ExtDir '/Anl_220202_1221_uniq_NormChr_Th5000'];

    PlotHists = 0; PlotFrags = 0;
%     HistIn.Genes = {'ACT1','CAT8','EFB1','GAL1','GAL2','GAL3','GAL4','GAL7','GAL80','GAS1','GCY1','HXT1','MAL31','MAL33','RPS14B','SUC2','TUB1'};
%     HistIn.Strns = {'WL93'}; HistIn.Conds = {'Dp1','Gp1'}; HistIn.Sorts = {'F','O'};
    HistIn.Genes = {'GAL1'}; % HistIn.Genes = {'GAL1','GAL2','GAL3','GAL7','GAL80'};
    HistIn.Strns = {'G'}; HistIn.Conds = {'Dp1','D0Gp1','Dm1G0','Dm2Gm1','Dm3Gm2','Dm4Gm3','Dm2Gm2','Dm2G0','Dm2Gp1','Gp1'}; HistIn.Sorts = {'O','F'};
    
        % choose only one of these two options:
    HistIn.ColorBySort = 0; HistIn.CondStyles = {'-','-','-','-','-','-','-','-'};
    HistIn.ColorByCond = 1; HistIn.SortStyles = {'-','--','.-'};
    
    PlotSums = 0; SumIn.Genes = HistIn.Genes; SumIn.Strns = HistIn.Strns; SumIn.Conds = HistIn.Conds; SumIn.Sorts = HistIn.Sorts;
    PlotCorels = 1;
    if strcmp(YFPExp,'L108')
        yfpOutDir = '../..\O25 - L108 Reanalysis\Stratedigm\L108 day1 day2 day4 day5\out_220202_2248'; %=with (220202_1044=wo) the added YFP data
            YFPIn.datadir = '../..\O25 - L108 Reanalysis\Stratedigm\L108 day1 day2 day4 day5\'; %SW25-26 YFP data was added from O19
            YFPIn.PlateMap = 'platemap_220202.xlsx'; YFPIn.muTH = 1.2;
            YFPIn.Plates = {'plate1','plate2','plate'}; % plate1 = rep1 (day1/4), plate2 = rep2 (day2/5), plate is from O19
            YFPIn.plot_hist=1; YFPIn.print_title=1; YFPIn.plot_means=1; YFPIn.plotgrad = 1; % choose gradient for PM plot 1=gluc/gal, 2=estradiol, 3=dox
%         CorelIn.Genes = {'GAL1'}; CorelIn.Strns = {'W','WL93'}; CorelIn.ConcatStrns = 1; % if concatStrns==1 make W the first in the lost
        CorelIn.Genes = {'GAL1'}; CorelIn.Strns = {'G','WL93'}; CorelIn.ConcatStrns = 0; % if concatStrns==1 make W the first in the lost
        CorelIn.Conds = {'Dp1','D0Gp1','Dm1G0','Dm2Gm1','Dm3Gm2','Dm4Gm3','Dm2Gm2','Dm2G0','Dm2Gp1','Gp1'}; CorelIn.Sorts = {'O','F'};
        CorelIn.ScaleATAC = 1; CorelIn.ErrBar = 1; CorelIn.YFPscale = 'log'; % set 'log' / 'linear' YScale
    elseif strcmp(YFPExp,'L135')
        yfpOutDir = '../../../..\GAL Decomposition\20190709 L135 DGs synGal4bsPr & GAL reporters\data\out_220201_1455';
            YFPIn.OldExtPath = '../../../..\GAL Decomposition\20190709 L135 DGs synGal4bsPr & GAL reporters\script\';
            YFPIn.datadir = '../../../..\GAL Decomposition\20190709 L135 DGs synGal4bsPr & GAL reporters\data\';
            YFPIn.PlateMap = 'platemap_220130_wo_GAL3_outlier.xlsx'; YFPIn.muTH = 1.2;
            YFPIn.plot_hist=0; YFPIn.print_title=1; YFPIn.plot_means=0; YFPIn.plotgrad = 1; % choose gradient for PM plot 1=gluc/gal, 2=estradiol, 3=dox
        CorelIn.Genes = {'GAL1','GAL2','GAL3','GAL7','GAL10','GAL80'}; CorelIn.Strns = {'W','WL93'}; % if concatStrns==1 make W the first in the lost
        CorelIn.ConcatStrns = 1; % concatenates strns into one (e.g. from diff Xlink experiments such as W and WL93)
        CorelIn.Conds = {'Dp1','D0Gp1','Dm1G0','Dm2Gm1','Dm3Gm2','Dm4Gm3','Gp1'}; CorelIn.Sorts = {'O','F'};
        CorelIn.ScaleATAC = 1; CorelIn.ErrBar = 1; CorelIn.YFPscale = 'linear'; % set 'log' / 'linear' YScale
    elseif strcmp(YFPExp,'L134') % (YFP from Ang's ALL GAL data)
        yfpOutDir = '../../../..\GAL Decomposition\20190708 L134 All GAL reporters from Ang\data\out_220203_1526';
            YFPIn.datadir = '../../../..\GAL Decomposition\20190708 L134 All GAL reporters from Ang\data\';
            YFPIn.PlateMap = 'platemap_220203.xlsx'; YFPIn.muTH = 1.2; YFPIn.DGformat = 1;
            YFPIn.Plates = {'plate1-yAL049-pGAL1','plate2-yAL139-pGAL2','plate3-yAL145-pGAL3','plate4-yAL142-pGAL7','plate5-yAL143-pGAL10','plate6-yAL144-pGAL80'};
            YFPIn.plot_hist=1; YFPIn.print_title=1; YFPIn.plot_means=0; YFPIn.plotgrad = 1; % choose gradient for PM plot 1=gluc/gal, 2=estradiol, 3=dox
        CorelIn.Genes = {'GAL1','GAL2','GAL3','GAL7','GAL10','GAL80'}; CorelIn.Strns = {'W','WL93'}; % if concatStrns==1 make W the first in the lost
        CorelIn.ConcatStrns = 1; % concatenates strns into one (e.g. from diff Xlink experiments such as W and WL93)
        CorelIn.Conds = {'Dp1','D0Gp1','Dm1G0','Dm2Gm1','Dm3Gm2','Dm4Gm3','Gp1'}; CorelIn.Sorts = {'O','F'};
        CorelIn.ScaleATAC = 1; CorelIn.ErrBar = 1; CorelIn.YFPscale = 'linear'; % set 'log' / 'linear' YScale
    end
end
%% O3 params
if strcmp(Exp,'O3')
%     datadir = '../..\O3 - 1st_ATACseq\data\align_200214_frag1000_match_refs/';
    datadir = '../..\O3 - 1st_ATACseq\data\align_220207_match_refs/';
    SmplNameON = 'S_WL93_Gp1_A'; % Used to calculate the scaled sum 'sumHistScaled'
    SmplNameOFF = 'S_WL93_Dp1_A'; % Used for both 'sumHistScaled' AND to plot histograms divided by the OFF histogram
    SmplNameONSNP = 'S_WS_Gp1_A';
    SmplNameOFFSNP = 'S_WS_Dp1_A'; 
    AnlIn.ReadThresh = 10000;
    AnlIn.BadSmpls = {};
    AnlIn.ut = 'uniq'; % 'uniq'/'tot' = remove/don't remove duplicates
    AnlIn.Norm = 'Loci'; % 'Loci'/'Chr'/'Uncap' = Normalize by control loci / chromosomes 3,11,14 / all uncaptured reads
    ExtDir = [datadir 'Ext_220207_2127_NF'];
        AnlDir = [ExtDir '/Anl_220207_2202_uniq_NormLoci_Th10000'];

    PlotHists = 0; PlotFrags = 0;
    HistIn.Genes = {'GAL1','HO_GAL1'};
    HistIn.Strns = {'WL93'}; HistIn.Conds = {'Dp1','Gp1'}; HistIn.Sorts = {'F','O'};
    % HistIn.Strns = {'GEV'}; HistIn.Sorts = {'A'}; HistIn.Conds = {'En','Ep1','Ep2','Ep3','Ep4','Ep5'}; % Estradiol SERIES O3
        % choose only one of these two options:
    HistIn.ColorBySort = 0; HistIn.CondStyles = {'-','-','-','-','-','-','-','-'};
    HistIn.ColorByCond = 1; HistIn.SortStyles = {'-','--','.-'};
    
    PlotSums = 0; SumIn.Genes = HistIn.Genes; SumIn.Strns = HistIn.Strns; SumIn.Conds = HistIn.Conds; SumIn.Sorts = HistIn.Sorts;
    PlotCorels = 1;
    % Corel HO_GAL1 (YFP and ATAC from O3)
    yfpOutDir = '../..\O3 - 1st_ATACseq\Stratedigm\20191017_GS-GSTT_W-WS_pure\data\out_220207_2233_w_bad_O3_samples';
        YFPIn.datadir = '../..\O3 - 1st_ATACseq\Stratedigm\20191017_GS-GSTT_W-WS_pure\data\'; %SW25-6 YFP data was added from O19
        YFPIn.PlateMap = 'platemap_w_bad_O3_samples.xlsx'; YFPIn.muTH = 1.2; % USE 'platemap.xlsx' TO ARTIFICIALLY REMOVE LOW READ ATAC SAMPLES
        YFPIn.Plates = {'plate1','plate2'}; % plate1 = W-WS pure sugars, plate2 = GS-GSTT gluc gradient
        YFPIn.plot_hist=1; YFPIn.print_title=1; YFPIn.plot_means=1; YFPIn.plotgrad = 1; % choose gradient for PM plot 1=gluc/gal, 2=estradiol, 3=dox
    CorelIn.Genes = {'HO_GAL1'}; CorelIn.Strns = {'GS','GSTT','WS'}; CorelIn.ConcatStrns = 0; % concatenates strns into one (e.g. from diff Xlink experiments such as W and WL93)
    CorelIn.Conds = {'Dp1','Gp1','D0Gm2','Dm1Gm2','Dm2Gm2','Dm3Gm2','Dm4Gm2','Dm5Gm2','Dm6Gm2','Gm2'}; CorelIn.Sorts = {'A'}; %,'U'
    CorelIn.ScaleATAC = 1; CorelIn.ErrBar = 1; CorelIn.YFPscale = 'linear'; % set 'log' / 'linear' YScale
end
%% O10 params
if strcmp(Exp,'O10')
    datadir = '../..\O10 - ATAC-seq GEV TATAdel\data\align_210205/';
    SmplNameON = 'S_WL93_Gp1_O_U'; % Used to calculate the scaled sum 'sumHistScaled'
    SmplNameOFF = 'S_WL93_Dp1_O_U'; % Used for both 'sumHistScaled' AND to plot histograms divided by the OFF histogram
    SmplNameONSNP = 'S_WS_Gp1_O_U';
    SmplNameOFFSNP = 'S_WS_Dp1_O_U';
    AnlIn.ReadThresh = 10000;
    AnlIn.BadSmpls = {};
    AnlIn.ut = 'uniq'; % 'uniq'/'tot' = remove/don't remove duplicates
    AnlIn.Norm = 'Loci'; % 'Loci'/'Chr'/'Uncap' = Normalize by control loci / chromosomes 3,11,14 / all uncaptured reads
    ExtDir = [datadir 'Ext_220207_2337_NF'];
        AnlDir = [ExtDir '/Anl_220207_2345_uniq_NormLoci_Th10000'];
    
    PlotHists = 0; PlotFrags = 0;
    HistIn.Genes = {'GAL1','HO_GAL1'};
    HistIn.Strns = {'WL93'}; HistIn.Conds = {'Dp1','Gp1'}; HistIn.Sorts = {'F','O'};
    % HistIn.Strns = {'GEV'}; HistIn.Sorts = {'O_A'}; HistIn.Conds = {'En','Ep1','Ep2','Ep3','Ep4','Ep5'}; % Estradiol SERIES O10
    % choose only one of these two options:
    HistIn.ColorBySort = 0; HistIn.CondStyles = {'-','-','-','-','-','-','-','-'};
    HistIn.ColorByCond = 1; HistIn.SortStyles = {'-','--','.-'};

    PlotSums = 0; SumIn.Genes = HistIn.Genes; SumIn.Strns = HistIn.Strns; SumIn.Conds = HistIn.Conds; SumIn.Sorts = HistIn.Sorts;
    PlotCorels = 1;
    % Corel HO_GAL1 (YFP and ATAC from O3)
    yfpOutDir = '../..\O3 - 1st_ATACseq\Stratedigm\20191017_GS-GSTT_W-WS_pure\data\out_220207_2233_w_bad_O3_samples';
        YFPIn.datadir = '../..\O3 - 1st_ATACseq\Stratedigm\20191017_GS-GSTT_W-WS_pure\data\'; %SW25-6 YFP data was added from O19
        YFPIn.PlateMap = 'platemap_w_bad_O3_samples.xlsx'; YFPIn.muTH = 1.2; % USE 'platemap.xlsx' TO ARTIFICIALLY REMOVE LOW READ ATAC SAMPLES
        YFPIn.Plates = {'plate1','plate2'}; % plate1 = W-WS pure sugars, plate2 = GS-GSTT gluc gradient
        YFPIn.plot_hist=1; YFPIn.print_title=1; YFPIn.plot_means=1; YFPIn.plotgrad = 1; % choose gradient for PM plot 1=gluc/gal, 2=estradiol, 3=dox
    CorelIn.Genes = {'HO_GAL1'}; CorelIn.Strns = {'GS','GSTT','WS'}; CorelIn.ConcatStrns = 0; % concatenates strns into one (e.g. from diff Xlink experiments such as W and WL93)
    CorelIn.Conds = {'Dp1','Gp1','D0Gm2','Dm1Gm2','Dm2Gm2','Dm3Gm2','Dm4Gm2','Dm5Gm2','Dm6Gm2','Gm2'}; CorelIn.Sorts = {'O_U'};
    CorelIn.ScaleATAC = 0; CorelIn.ErrBar = 0; CorelIn.YFPscale = 'linear'; % set 'log' / 'linear' YScale
end
%% O16 params
if strcmp(Exp,'O16')
    % datadir = '../..\O16 - ATACseq O11+O15\data\align_210526/';
    datadir = '../..\O16 - ATACseq O11+O15\data\align_210526_del/'; % I artificially deleted the GAL1_HO reads for samples 150, 188, 192 that were mistakenly aligned to YeastGenome_HO_TTDEL
    SmplNameON = 'S_WL93_Gp1_H_m'; % Used to calculate the scaled sum 'sumHistScaled'
    SmplNameOFF = 'S_WL93_Dp1_H_m'; % Used for both 'sumHistScaled' AND to plot histograms divided by the OFF histogram
    SmplNameONSNP = SmplNameON; %workaround cause there are no SNPs in this run
    SmplNameOFFSNP = SmplNameOFF;
    AnlIn.ReadThresh = 10000;
    AnlIn.BadSmpls = {};
    AnlIn.ut = 'uniq'; % 'uniq'/'tot' = remove/don't remove duplicates
    AnlIn.Norm = 'Loci'; % 'Loci'/'Chr'/'Uncap' = Normalize by control loci / chromosomes 3,11,14 / all uncaptured reads
    ExtDir = [datadir 'Ext_220205_2240_NF'];
        AnlDir = [ExtDir '/Anl_220205_2243_uniq_Loci_Th10k'];
        
    PlotHists = 0; PlotFrags = 0;
    HistIn.Genes = {'GAL1'}; % HistIn.Genes = {'GAL1','GAL2','GAL3','GAL7','GAL80'};
    HistIn.Strns = {'G'}; HistIn.Conds = {'Dp1','D0','D0Gm2','Dm1Gm2','Dm2Gm2','Dm3Gm2','Dm4Gm2','Gm2','Gp1'}; HistIn.Sorts = {'H_m'};
        % choose only one of these two options:
    HistIn.ColorBySort = 0; HistIn.CondStyles = {'-','-','-','-','-','-','-','-'};
    HistIn.ColorByCond = 1; HistIn.SortStyles = {'-','--','.-'};
    
    PlotSums = 0; SumIn.Genes = HistIn.Genes; SumIn.Strns = HistIn.Strns; SumIn.Conds = HistIn.Conds; SumIn.Sorts = HistIn.Sorts;
    PlotCorels = 1; PLATE = 'Plate1';
    if strcmp(PLATE,'Plate2')
        yfpOutDir = '../..\O11 - Xlink chromatin mutants\O11cd - 20210314 BSSG Xlink gal80del bg mutants at later times\data\out_220206_2144_wo_outlier';
            YFPIn.PlateMap = 'platemap_220206_wo_outlier.xlsx';
            YFPIn.datadir = '../..\O11 - Xlink chromatin mutants\O11cd - 20210314 BSSG Xlink gal80del bg `mutants at later times\data\';
            YFPIn.Plates = {'t45D','t68D','t110C','plate'};
            YFPIn.plot_hist=1; YFPIn.print_title=1; YFPIn.plot_means=1; YFPIn.plotgrad = 1; % choose gradient for PM plot 1=gluc/gal, 2=estradiol, 3=dox
        CorelIn.Genes = {'GAL1'}; CorelIn.Strns = {'G','GSNF','GPGD','GSPT','WL93'};
        CorelIn.ConcatStrns = 0; % if concatStrns==1 make W the first in the lost
        CorelIn.Conds = {'Dp1','D0','D0Gm2','Dm1Gm2','Dm2Gm2','Dm3Gm2','Dm4Gm2','Gm2','Gp1'};
        CorelIn.TechReps = 1; CorelIn.Sorts = {'H_m'};
        CorelIn.ScaleATAC = 0; CorelIn.ErrBar = 0; CorelIn.YFPscale = 'log'; % set 'log' / 'linear' YScale
    elseif strcmp(PLATE,'Plate1')
        yfpOutDir = '../..\O15 - Xlink GEV gradient + WT pure gal-gluc\O15b Xlink GEV + WT\data\out_220222_1045';
            YFPIn.PlateMap = 'platemap_220222.xlsx';
            YFPIn.datadir = '../..\O15 - Xlink GEV gradient + WT pure gal-gluc\O15b Xlink GEV + WT\data\';
            YFPIn.Plates = {'t8','plate'};
            YFPIn.plot_hist=1; YFPIn.print_title=1; YFPIn.plot_means=1; YFPIn.plotgrad = 2; % choose gradient for PM plot 1=gluc/gal, 2=estradiol, 3=dox
        CorelIn.Genes = {'GAL1'}; CorelIn.Strns = {'WL93','GEVDD','GEV'}; CorelIn.ConcatStrns = 0; % if concatStrns==1 make W the first in the lost
        CorelIn.Conds = {'Gp1','Dp1','Dp1Ep01','Dp1Ep02','Dp1Ep03','Dp1Ep04','Dp1Ep05','Dp1Ep06','Dp1Ep07','Dp1Ep08','Dp1Ep09'};
        CorelIn.TechReps = 1; CorelIn.Sorts = {'H_U'};
        CorelIn.ScaleATAC = 1; CorelIn.ErrBar = 1; CorelIn.YFPscale = 'log'; % set 'log' / 'linear' YScale
    end
end
%% O19 params
if strcmp(Exp,'O19')
    datadir = '../..\O19 - Test SPRI_cell_Tn5 ratios\data\align_211110/';
    SmplNameON = 'S_W_Gp1_HH_U'; % Used to calculate the scaled sum 'sumHistScaled'
    SmplNameOFF = 'S_W_Dp1_HH_U'; % Used for both 'sumHistScaled' AND to plot histograms divided by the OFF histogram
    SmplNameONSNP = SmplNameON; %workaround cause there are no SNPs in this run
    SmplNameOFFSNP = SmplNameOFF;
    AnlIn.BadSmpls = {};
    % AnlIn.BadSmpls = {AnlIn.BadSmpls 'O19_S_W_Dp1_HH_M1_2'}; % S104 - outlier seen in fragment size distribution
    % AnlIn.BadSmpls = {AnlIn.BadSmpls 'O19_S_W_Gp1_HH_M3_2'}; % S59 - outlier seen in Tn5-insertion histogram
    % AnlIn.BadSmpls = {AnlIn.BadSmpls 'O19_S_W_Dp1_HH_MM1_1'}; % S69 - outlier seen in Tn5-insertion histogram
    AnlIn.ReadThresh = 10000;
    AnlIn.ut = 'uniq'; % 'uniq'/'tot' = remove/don't remove duplicates
    AnlIn.Norm = 'Chr'; % 'Loci'/'Chr'/'Uncap' = Normalize by control loci / chromosomes 3,11,14 / all uncaptured reads
    ExtDir = [datadir 'Ext_211111_0855_NF'];
        AnlDir = [ExtDir '/Anl_211111_1339_tot_Th10000'];
    
    PlotHists = 0; PlotFrags = 0;
    HistIn.Genes = {'GAL1'};
    HistIn.Strns = {'W'}; HistIn.Conds = {'Dp1','Gp1'}; HistIn.Sorts = {'HH_U'};
        % choose only one of these two options:
    HistIn.ColorBySort = 0; HistIn.CondStyles = {'-','-','-','-','-','-','-','-'};
    HistIn.ColorByCond = 1; HistIn.SortStyles = {'-','--','.-'};
    
    PlotSums = 0; SumIn.Genes = HistIn.Genes; SumIn.Strns = HistIn.Strns; SumIn.Conds = HistIn.Conds; SumIn.Sorts = HistIn.Sorts;
    PlotCorels = 1;
%     yfpOutDir = '../..\O19 - Test SPRI_cell_Tn5 ratios\210721 O19a Xlink\out_220202_1044';
        YFPIn.datadir = '../..\O19 - Test SPRI_cell_Tn5 ratios\210721 O19a Xlink\';
        YFPIn.PlateMap = 'platemap_220202.xlsx';
        YFPIn.Plates = {'plate'};
        YFPIn.plot_hist=1; YFPIn.print_title=1; YFPIn.plot_means=1; YFPIn.plotgrad = 1; % choose gradient for PM plot 1=gluc/gal, 2=estradiol, 3=dox
    CorelIn.Genes = {'GAL1'}; CorelIn.Strns = {'W'}; CorelIn.ConcatStrns = 0; % if concatStrns==1 make W the first in the lost
    CorelIn.Conds = {'Dp1','Gp1'}; CorelIn.Sorts = {'HH_U'};
    CorelIn.YFPscale = 'log'; % set 'log' / 'linear' YScale
end
%% O24 params
if strcmp(Exp,'O24')
    datadir = '../..\O24 - ATAC-seq O14f\data\align_220208/'; % I artificially deleted the GAL1_HO reads for samples 150, 188, 192 that were mistakenly aligned to YeastGenome_HO_TTDEL
    SmplNameON = 'S_WO19_Gp1_AM_P'; % Used to calculate the scaled sum 'sumHistScaled'
    SmplNameOFF = 'S_WO19_Dp1_AM_P'; % Used for both 'sumHistScaled' AND to plot histograms divided by the OFF histogram
    %%%%%%% ALSO TRY THE 3 WL93_Gp1 SAMPLES AS NORMALIZING (there are 5 WO19_Gp1 samples) %%%%%%%%%
    SmplNameONSNP = SmplNameON; %workaround cause there are no SNPs in this run
    SmplNameOFFSNP = SmplNameOFF;
    AnlIn.ReadThresh = 10000;
%     AnlIn.BadSmpls = {'S_ZF1MED6_X320_AM_P3','S_ZF1TTI1_X0_AM_P1','S_ZF2TTI1_X320_AM_P1_2'};
    AnlIn.BadSmpls = {'S_ZF1MED6_X320_AM_P3','S_ZF1TTI1_X0_AM_P1'};
    AnlIn.ut = 'uniq'; % 'uniq'/'tot' = remove/don't remove duplicates
    AnlIn.Norm = 'Loci'; % 'Loci'/'Chr'/'Uncap' = Normalize by control loci / chromosomes 3,11,14 / all uncaptured reads
    ExtDir = [datadir 'Ext_220209_1640_NF_no_SPRI']; %wo SPRI information
        AnlDir = [ExtDir '/Anl_220209_2156_uniq_Loci_Th10k'];
        
    PlotHists = 0; PlotFrags = 0;
%     HistIn.Genes = {'ZF1_CYC1'};
%     HistIn.Strns = {'ZF1','ZF1TTI1','ZF1HDA3','ZF1SWD3','ZF1MED6'};
    HistIn.Genes = {'ZF2_CYC1'};
    HistIn.Strns = {'ZF2','ZF2HDA3','ZF2MED6','ZF2TTI1','ZF2SWD3'};
    HistIn.Conds = {'X0','X40','X80','X160','X320','X640','X1280','X2560'}; HistIn.Sorts = {'AM_P','HB_P'};
        % choose only one of these two options:
    HistIn.ColorBySort = 0; HistIn.CondStyles = {'-','-','-','-','-','-','-','-'};
    HistIn.ColorByCond = 1; HistIn.SortStyles = {'-','--','.-'};
    
    PlotSums = 0; SumIn.Genes = HistIn.Genes; SumIn.Strns = HistIn.Strns; SumIn.Conds = HistIn.Conds; SumIn.Sorts = HistIn.Sorts;
    SumIn.TechReps = 1:2;

    PlotCorels = 1;
    if strcmp(YFPExp,'O24')
        yfpOutDir = '../..\O14 - Xlink Khalil library\O14f Xlink 8 strains + WT\data\out_220212_2228';
            YFPIn.PlateMap = 'platemap_220211.xlsx';
            YFPIn.datadir = '../..\O14 - Xlink Khalil library\O14f Xlink 8 strains + WT\data\';
            YFPIn.Plates = {'P1','P2','P3','plateL93','plateO19'}; %YFPIn.Plates = {'P2','P3'}; 
            YFPIn.plot_hist=1; YFPIn.print_title=1; YFPIn.plot_means=1; YFPIn.plotgrad = 3; % choose gradient for PM plot 1=gluc/gal, 2=estradiol, 3=dox
        CorelIn.Genes = {'ZF2_CYC1'}; CorelIn.Strns = {'ZF2MED6','ZF2SWD3'}; %,'ZF2HDA3','ZF2TTI1'
        CorelIn.ConcatStrns = 0; % if concatStrns==1 make W the first in the lost
%         CorelIn.TechReps = 1:2; CorelIn.Conds = {'X0','X40','X80','X160','X320','X640','X1280','X2560'}; CorelIn.Sorts = {'HB_P','AM_P'}; %
        CorelIn.TechReps = 1:2; CorelIn.Conds = {'X0','X40','X80','X160'}; CorelIn.Sorts = {'L_P'}; %,'X320','X640','X1280','X2560'
        CorelIn.ScaleATAC = 0; CorelIn.ErrBar = 1; CorelIn.YFPscale = 'linear'; % set 'log' / 'linear' YScale / any other entry would not fit at all
    end
end
%% Define summing regions and binding sites
% I defined the below regions by looking for inaccessible regions in 2% gluc with 3bp smoothing window (but should shift ends by 5bp to make sure)
Regs.GAL1 = struct('Nuc_m2',[278390,278520],'NFR',[278521,278729],'NUC_NFR',[278575,278677],'Nuc_m1',[278730,278860],'Nuc_p1',[278917,279020],'Prom',[278353,279020]);
Regs.HO_GAL1 = struct(                      'NFR',[3216,3378],    'NUC_NFR',[3260,3360],'Nuc_m1',[3075,3205],    'Nuc_p1',[2908,3018],    'Prom',[2911,3378]);
Regs.GAL2 = struct('Nuc_m1',[289930,290050],'Nuc_p1',[290101,290211]);
Regs.GAL3 = struct('Nuc_m2',[463000,463100],'Nuc_m1',[463193,463279]); % Nuc+1 doesn't show difference between gal and gluc
Regs.GAL7 = struct('Nuc_m1',[275827,275934],'Nuc_p1',[275530,275660]);
Regs.GAL10 = struct('Nuc_m2',[278730,278860],'Nuc_m1',[278390,278510],'NFR',[278521,278729],'NUC_NFR',[278575,278677]);
Regs.GAL80 = struct('Nuc_p1',[171522,171621]); % old: 'Nuc_p1',[171490,171594]
Regs.ZF1_CYC1 = struct('Nuc_m1',[1658,1779],'Prom',[826,1929]);
Regs.ZF2_CYC1 = struct('Nuc_m1',[1658,1779],'Prom',[826,1953]);

% Centers from Chereji et al 2018:
% GAL1_Nuc-1 = 278471, GAL1_NDR = 278631 (width 172), GAL1_Nuc+1 = 278790, GAL1_TSS = 279021, GAL1_TTS = 280655
% GAL10_TSS = 278352, GAL10_TTS = 276151

% Nuc centers from V plot:
% GAL1_Nuc-1 = 278789, GAL1_Nuc+1 = 278930
% GAL2_Nuc-1 = 289985 / 290004 / 289993, GAL2_Nuc+1 = 290154
% GAL7_Nuc+1 = 275587
% GAL80_Nuc+1 = 171578
% GAL10_Nuc+1 = 278447

Sites.Gal4bs.GAL1 = {[278568,278584],[278587,278603],[278605,278621],[278667,278691]}; %[278567,278586,278604,278668];
Sites.Mig1bs.GAL1 = {[278751,278763],[278821,278832]}; %[278750,278820];
Sites.TATA.GAL1 = [278874,278881];
Sites.TATA.HO_GAL1 = [3053,3060];
Sites.TSS.GAL1 = [278955,278960];
Sites.ZFbs.ZF1_CYC1 = {[1663,1673]};
Sites.ZFbs.ZF2_CYC1 = {[1663,1673],[1688,1698]};
Sites.TATA.ZF1_CYC1 = [1747,1809]; % I COMBINED THESE TWO REGIONS INTO ONE [1747,1755],[1801,1809]
Sites.TATA.ZF2_CYC1 = [1772,1834]; % I COMBINED THESE TWO REGIONS INTO ONE [1772,1780],[1826,1834]

% OFFSET=commonStart(1);
StrnColors = struct('S_W',[0 0 0],'S_WL93',[0.5 0.5 0.5],'S_WS',[0 0 1],'S_GS',[0.64 0.08 0.18],'S_GSTT',[0.4660 0.6740 0.1880],...
    'S_M',[1 0 0],'S_MG',[0 1 0],'S_TETM',[0 0 1],'S_GEV',[0 0.5 0],'S_GEVDD',[0.75 0 0],...
    'S_G',[0 0 1],'S_GSNF',[0.87 0.49 0],'S_GPGD',[0 0 0],'S_GSPT',[0 1 1],...
    'S_ZF1TTI1',[0 0 1],'S_ZF1HDA3',[0.87 .49 0],'S_ZF1SWD3',[0 0 0],'S_ZF1MED6',[0 1 0],'S_ZF1',[0.5 0.5 0.5],...
    'S_ZF2TTI1',[0 0 1],'S_ZF2HDA3',[0.87 .49 0],'S_ZF2SWD3',[0 0 0],'S_ZF2MED6',[0 1 0],'S_ZF2',[0.5 0.5 0.5]);
RepMarkers = {'d','o','s','p'}; % diamond, circle, square?, pentagon? for bioreps 1,2,3,4
%% Extract
% Read ATACseq txt files, based on IDglobal.xlsx/Samples and turn it to structs

if exist('ExtDir','var')
    load([ExtDir '/ExtOut.mat'],'ExtOut') %Load the output from the extraction script
else
    ExtIn.MAPQmin = 10;
    [ExtDir,ExtOut] = a_Extract_ATAC_txtfiles(ExtIn,Exp,datadir);
    
    CurrentCodeFile = [mfilename '.m'];
    [status,msg] = copyfile(CurrentCodeFile,ExtDir);
    if status==0
        error('Error. pipeline code file was not copied after extraction step')
    end
end
%% Arrange
% Arrange insertions histograms and their sums based on IDglobal.xlsx/Table
% into a strcuct of gene X strain X condition X sort gate

if exist('AnlDir','var')
    load([AnlDir '/AnlOut.mat'],'AnlOut') %Load the output from the analysis script
else
    AnlIn.checkSmplCount = 1; %makes sure that SmplCount in AnlOut and ExtOut are the same for each gene
    
    if strcmp(AnlIn.Norm,'Loci')
        if strcmp(AnlIn.ut,'uniq')
            load([ExtDir '/Hists.mat'],'HistUniqNormLoci') %Load hists from the extraction script
            Hist = HistUniqNormLoci; clear HistUniqNormLoci
        elseif strcmp(AnlIn.ut,'tot')
            load([ExtDir '/Hists.mat'],'HistTotNormLoci') %Load hists from the extraction script
            Hist = HistTotNormLoci; clear HistTotNormLoci
        end
    elseif strcmp(AnlIn.Norm,'Chr')
        if strcmp(AnlIn.ut,'uniq')
            load([ExtDir '/Hists.mat'],'HistUniqNormChr') %Load hists from the extraction script
            Hist = HistUniqNormChr; clear HistUniqNormChr
        elseif strcmp(AnlIn.ut,'tot')
            load([ExtDir '/Hists.mat'],'HistTotNormChr') %Load hists from the extraction script
            Hist = HistTotNormChr; clear HistTotNormChr
        end
    end
    
    [AnlDir,AnlOut] = a_Analyze_ATAC(Exp,datadir,ExtDir,ExtOut,Hist,AnlIn,SmplNameON,SmplNameOFF,SmplNameONSNP,SmplNameOFFSNP,Regs);
    
    CurrentCodeFile = [mfilename '.m'];
    [status,msg] = copyfile(CurrentCodeFile,AnlDir);
    if status==0
        error('Error. pipeline code file was not copied after analysis step')
    end
end
%% Plot Hists for chosen samples       %%%%%%%%%% MOVE EXP SPECIFIC PARAMS TO THE TOP %%%%%%%%%%%%%%%%
HistIn.Legend = [0 0 1 0 0 0]; %Collect these for legend [TYPE,STRN,COND,TAGMENT,CELLNUM,EXP]
HistIn.LegendOnOff = [0 1 1 0 0 0]; %Collect these for ON/OFF smpl legend [TYPE,STRN,COND,TAGMENT,CELLNUM,EXP]
HistIn.DirName = ['Hist_' HistIn.Strns{:} '_' HistIn.Sorts{:}];
if length(HistIn.DirName)>20
    HistIn.DirName = 'Hist_TOOLONG';
end
HistIn.Colors = {[0.6350 0.0780 0.1840],[0.4660 0.6740 0.1880],[0 0 0],[0.8500 0.3250 0.0980],[0 0.4470 0.7410],[0.4940 0.1840 0.5560],[0.3010 0.7450 0.9330],[0.9290 0.6940 0.1250],[0.6350 0.0780 0.1840],[0.4660 0.6740 0.1880],[0 0 0],[0.8500 0.3250 0.0980]};
HistIn.plotOnOffHists = 0; HistIn.plotlegend = 1;
HistIn.AvgRep = 1; HistIn.AvgTechRep = 0; HistIn.smoothWidth = 10; HistIn.savefigs = 1;

HistIn.DivByHist = 0; %HistIn.HistDivStrn = 'S_WL93'; HistIn.HistDivCond = 'Gp1'; HistIn.HistDivSort = 'H_U';

if PlotHists==1
    if HistIn.savefigs==1
        HistIn.Dir = [AnlDir '/' HistIn.DirName '_Avg' num2str(HistIn.AvgRep) '_Smooth' num2str(HistIn.smoothWidth) '_' datestr(now,'yymmdd_HHMM')];% to save the hists
        mkdir(HistIn.Dir);
        
        CurrentCodeFile = [mfilename '.m'];
        [status,msg] = copyfile(CurrentCodeFile,HistIn.Dir);
        if status==0
            error('Error. pipeline code file was not copied after histogram plotting step')
        end
    end
    a_plot_ATAC_Hists(Exp,datadir,HistIn,AnlIn,AnlOut,Regs,Sites);
end

%% Plot fragment lengths and length vs mid point for chosen samples
FragIn.Genes = HistIn.Genes; FragIn.Strns = HistIn.Strns; FragIn.Conds = HistIn.Conds; FragIn.Sorts = HistIn.Sorts;
% FragIn.Conds = {'Dp1'}; FragIn.Sorts = {'HH_M','H_M','L_M','LL_M'}; %FragIn.Reps = [1]; %[2 4 6];

FragIn.savefigs = 1; FragIn.stackStrns = 0;
FragIn.plotlenVSmid = 0; FragIn.lenVSmidBioRep = 1; FragIn.lenVSmidTechRep = 1;%plot only fragment length vs midpoint of chosen rep
FragIn.Colors = HistIn.Colors; FragIn.CondStyles = HistIn.CondStyles;
FragIn.ColorBySort = 1; FragIn.plotlegend = 1; FragIn.AvgRep = 0; FragIn.AvgTechRep = 0;

if FragIn.plotlenVSmid == 1
    FragIn.DirName = ['lenVSmid_' FragIn.Genes{:} '_' FragIn.Strns{:} '_' FragIn.Conds{:} '_' FragIn.Sorts{:}];
else
    FragIn.DirName = ['Frag_' FragIn.Genes{:} '_' FragIn.Strns{:} '_' FragIn.Conds{:} '_' FragIn.Sorts{:}];
end
if length(FragIn.DirName) > 20
    FragIn.DirName = 'Frag_TOOLONG';
end

if PlotFrags == 1
    if FragIn.savefigs == 1
        FragIn.Dir = [AnlDir '/' FragIn.DirName '_' datestr(now,'yymmdd_HHMM')];% to save the figs
        mkdir(FragIn.Dir);
        
        CurrentCodeFile = [mfilename '.m'];
        [status,msg] = copyfile(CurrentCodeFile,FragIn.Dir);
        if status==0
            error('Error. pipeline code file was not copied after histogram plotting step')
        end
    end
    a_plot_ATAC_Frags(Exp,FragIn,AnlOut,Regs,Sites);
end

%% plot sums of chosen samples and all regions   %%%%%%%%%% MOVE EXP SPECIFIC PARAMS TO THE TOP %%%%%%%%%%%%%%%%
if PlotSums==1
    SumIn.Xlabel = [0 0 1 1 0 0]; %Collect these for legend [TYPE,STRN,COND,TAGMENT,CELLNUM,EXP]
    SumIn.DirName = ['Sum_' SumIn.Genes{:} '_' SumIn.Strns{:} '_' SumIn.Sorts{:}];
    if length(SumIn.DirName)>20
        SumIn.DirName = 'Sum_TOOLONG';
    end
    SumIn.Colors = StrnColors;
    SumIn.Markers = RepMarkers;
    SumIn.savefigs = 1; SumIn.Norm = 0;
    
    if SumIn.savefigs==1
        SumIn.Dir = [AnlDir '/' SumIn.DirName '_' datestr(now,'yymmdd_HHMM')];% to save the Sums
        mkdir(SumIn.Dir);
        
        CurrentCodeFile = [mfilename '.m'];
        [status,msg] = copyfile(CurrentCodeFile,SumIn.Dir);
        if status==0
            error('Error. pipeline code file was not copied after sum plotting step')
        end
    end
    
    SumOut = a_plot_ATAC_Sums(Exp,datadir,SumIn,AnlOut,Regs);
end

%% extract YFP and plot correlations
if PlotCorels==1
    %% Extract expression data      %%%%%%%%% SEE TOP OF CODE FOR EXP SPECIFIC PARAMS %%%%%%%%
    YFPIn.Colors = StrnColors;
    YFPIn.Markers = RepMarkers;
    YFPIn.saveoutput=1; YFPIn.filDeb=0; %YFPIn.DebGates = 'debris manualgate plate1 E05 20190716.mat';
    if exist('yfpOutDir','var')
        load([yfpOutDir '/extracted metrics.mat'],'logYFP') %Load the output from the analysis script
        if exist('yfpOutDir2','var')
            logYFP1 = logYFP;
            load([yfpOutDir2 '/extracted metrics.mat'],'logYFP') %Load the output from the analysis script
            for iGene = fieldnames(logYFP1)'
                for iStrn = fieldnames(logYFP1.(iGene{1}))'
                    for iCond = fieldnames(logYFP1.(iGene{1}).(iStrn{1}))'
                        % assign logYFP1 into logYFP ASSUMING they have distinct (gene X strn X cond) sets
                        logYFP.(iGene{1}).(iStrn{1}).(iCond{1}) = logYFP1.(iGene{1}).(iStrn{1}).(iCond{1});
                    end
                end
            end
        end
    else
        if isfield(YFPIn,'OldExtPath')
            logYFP = a_extract_YFP_from_old_metrics(YFPIn);
        else
            logYFP = a_extract_YFP_pm_IndWells(YFPIn);
        end
    end
    
    %% Plot correlations          %%%%%%%%% SEE TOP OF CODE FOR EXP SPECIFIC PARAMS %%%%%%%%
    CorelIn.Colors = StrnColors;
    CorelIn.Markers = RepMarkers;
    CorelIn.DirName = ['Corel_' YFPExp '_' CorelIn.Strns{:}];
    if length(CorelIn.DirName)>25
        CorelIn.DirName = ['Corel_' YFPExp '_TOOLONG'];
    end
    CorelIn.savefigs = 1;
    
    if CorelIn.savefigs==1
        if CorelIn.ScaleATAC == 0
            if CorelIn.ErrBar == 0
                CorelIn.Dir = [AnlDir '/' CorelIn.DirName '_' datestr(now,'yymmdd_HHMM') '_' CorelIn.YFPscale(1:3)];
            elseif CorelIn.ErrBar == 1
                CorelIn.Dir = [AnlDir '/' CorelIn.DirName '_' datestr(now,'yymmdd_HHMM') '_' CorelIn.YFPscale(1:3) '_err'];
            end
        elseif CorelIn.ScaleATAC == 1
            if CorelIn.ErrBar == 0
                CorelIn.Dir = [AnlDir '/' CorelIn.DirName '_' datestr(now,'yymmdd_HHMM') '_' CorelIn.YFPscale(1:3) '_scl'];
            elseif CorelIn.ErrBar == 1
                CorelIn.Dir = [AnlDir '/' CorelIn.DirName '_' datestr(now,'yymmdd_HHMM') '_' CorelIn.YFPscale(1:3) '_scl_err'];
            end
        end
        mkdir(CorelIn.Dir);
        CurrentCodeFile = [mfilename '.m'];
        [status,msg] = copyfile(CurrentCodeFile,CorelIn.Dir);
        if status==0
            error('Error. pipeline code file was not copied after Correlation plotting step')
        end
    end
    a_plot_ATAC_YFP_corel(Exp,YFPExp,CorelIn,logYFP,AnlOut,Regs);
end