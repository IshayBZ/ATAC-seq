function [ExtDir,ExtOut] = Extract_ATAC_txtfiles(ExtIn,Exp,datadir)
% PURPOSE: extract structures of ATACseq start site histograms, positions using the functions 'readTableIn.m',
%     'uploadAtacSeqData.m', and 'makeStartHistograms3.m'. Outputs the file 'tables_and_histograms.mat'
% NOTE: expects a folder 'data' that contains experiment-specific subfolders with .txt files and 
%     IDglobal.csv file that contains:
    % a. sample # to sample name conversion for that experiment
    % b. coordinates for captured genes
    % c. table sheet of strains and conditions in experiment
    % d. key sheet that converts short names to descriptions

%% Extract the atacHist, atacStartStop structs (collected in umbrella structs)
tic

[UniqCaptLR,UniqAllLR,TotCaptLR,TotAllLR,BadReads,ExtOut.SmplCount] = uploadATACseqData(Exp,datadir,ExtIn.MAPQmin);

ExtOut.outTab = readTableIn(UniqCaptLR,UniqAllLR,TotCaptLR,TotAllLR,BadReads);

[HistUniqNormLoci,HistUniqNormChr,HistTotNormLoci,HistTotNormChr,ExtOut.StartStop,ExtOut.FragLenHist] ...
    = makeTn5CutsHists(UniqCaptLR,TotCaptLR,ExtOut.outTab);
toc
%160sec - optimize speed if desired, but ok for now (200110)

%% SAVE OUTPUT IN A NEW FOLDER, TOGETHER WITH A SNAPSHOT OF THE CODE
ExtDir = [datadir 'Ext_' datestr(now,'yymmdd_HHMM') '_NF']; % NF = new format (was made on 211028)
mkdir(ExtDir);
save([ExtDir '/Reads.mat'],'UniqCaptLR','UniqAllLR','TotCaptLR','TotAllLR','BadReads');
save([ExtDir '/Hists.mat'],'HistUniqNormLoci','HistUniqNormChr','HistTotNormLoci','HistTotNormChr');
save([ExtDir '/ExtOut.mat'],'ExtOut','Exp');
CurrentCodeFile = [mfilename '.m'];
copyfile(CurrentCodeFile,ExtDir);