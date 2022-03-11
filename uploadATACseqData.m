function [UniqCaptLR,UniqAllLR,TotCaptLR,TotAllLR,BadReads,SmplCount] = uploadATACseqData(Exp,datadir,MAPQmin)
% grab start and stop positions for correct alignments, using samfile columns written to txt files
% CaptData(:,1) is samfile $4 (POS) - if TLEN>0 this is the left end of a paired-end fragment
% CaptData(:,2) is now length($10) which means nothing
% CaptData(:,3) is $5 (MAPQ) - the mapping quality, e.g. MAPQ>10 means Prob(correct mapping)>0.9
% CaptData(:,4) is $9 (TLEN) - the signed distance between lefmost and rightmost mapped bases (+ for F, - for R, 0 for unpaired read)
% (i.e. if TLEN>0 then POS+TLEN is the right end of a paired-end fragment + 1)
% (NOTE: if uniq==1 we're throwing unpaired reads because we can't tell if they're PCR duplicates)

%This function expects one directory in which you have put all the files
% UPDATED 20200127 - to correctly read the right end of fragment from the txt file:
                     % - to use column 1+4 instead of 2 to get right end of fragment
                     % - to use only MAPQ>10 (i.e. Prob(correct mapping)>0.9)
                     % - to allow unpaired reads as well for the non-unique data
% UPDATED 20210209 - to also extract the mitochondrial chromosome and HO locus
% UPDATED 20210407 - to extract fragment length and midpoint
% UPDATED 20211027 - to save only LR (left & right ends) to minimize file sizes (farment length and midpoint will be calculated from that)

%load Look up table - This expects a file called ['IDglobal_' 'Exp' '.xlsx'] in the dirPath.
try
    num2strain = readtable([datadir '../../IDglobal_',Exp,'.xlsx'],'Sheet','Samples');
catch
    num2strain = readtable([datadir 'IDglobal_',Exp,'.xlsx'],'Sheet','Samples');
end
SmplCount = size(num2strain,1);

%Grabs txt of captured data
files = dir([datadir,filesep,'*.txt']);
filenames = {files.name};
filenames = natsort(filenames);
for iFile=1:numel(filenames)
    curName = filenames(iFile);
    nameTokens = regexp(curName{1},'starts_','split');
    geneName = nameTokens{1};
    sampleNum = str2num(nameTokens{2}(2:strfind(nameTokens{2},'.txt')-1));
    sampleName = num2strain{num2strain{:,1}==sampleNum,2};
    CaptData = importdata([datadir,filesep,curName{1}]); %Load Data
    
    %% Continue with Tn5CutsT and U for tot and uniq reads
    if isempty(CaptData) %for samples with no reads
        BadReads.Capt.LowMapQ.(geneName).(sampleName{1}) = [];
        BadReads.Capt.Unpair.(geneName).(sampleName{1}) = [];
        BadReads.Capt.FarPairs.(geneName).(sampleName{1}) = [];
        TotCaptLR.(geneName).(sampleName{1}) = [];
        UniqCaptLR.(geneName).(sampleName{1}) = [];        
        continue
    end
    
    %saving relevant half of reads which didn't map well enough
    BadMapCaptData = CaptData(CaptData(:,3)<=MAPQmin & CaptData(:,4)>0,1:4);
    BadReads.Capt.LowMapQ.(geneName).(sampleName{1}) = [BadMapCaptData(:,1) BadMapCaptData(:,1) + BadMapCaptData(:,4)];
    CaptData = CaptData(CaptData(:,3)>MAPQmin,1:4); %keeping reads which mapped well
    
    % saving unpaired reads and pairs that are too far
    BadReads.Capt.Unpair.(geneName).(sampleName{1}) = CaptData(CaptData(:,4)==0,1);
    BadReads.Capt.FarPairs.(geneName).(sampleName{1}) = [CaptData(CaptData(:,4)>20000,1) CaptData(CaptData(:,4)>20000,1) + CaptData(CaptData(:,4)>20000,4)];
    
    %keep paired reads where col 1 is left (col4>0) and right is not too far
    CaptData = CaptData(CaptData(:,4)>0 & CaptData(:,4)<20000,1:4);
    
    % Saving left & right ends (of paired reads)
    TotCaptLR.(geneName).(sampleName{1}) = [CaptData(:,1) CaptData(:,1) + CaptData(:,4)];
    
    % Saving unique left & right ends (of paired reads)
    [~,ia,~] = unique(TotCaptLR.(geneName).(sampleName{1}),'rows');
    UniqCaptData = CaptData(ia,:);
    UniqCaptLR.(geneName).(sampleName{1}) = [UniqCaptData(:,1) UniqCaptData(:,1) + UniqCaptData(:,4)];
end

%% Extracting data from all chromosomes
%Grabs txt of All chr data
datadirAll = [datadir 'AllChr'];
filesAll = dir([datadirAll,filesep,'*.txt']);
filenamesAll = {filesAll.name};
filenamesAll = natsort(filenamesAll);
for iFile=1:numel(filenamesAll)
    curName = filenamesAll(iFile);
    nameTokens = regexp(curName{1},'starts_','split');
    chrName = nameTokens{1};
    sampleNum = str2num(nameTokens{2}(2:strfind(nameTokens{2},'.txt')-1));
    sampleName = num2strain{num2strain{:,1}==sampleNum,2};
    AllData = importdata([datadirAll,filesep,curName{1}]); %Load Data
    
    if isempty(AllData) % for samples with no reads
        BadReads.All.LowMapQ.(chrName).(sampleName{1}) = [];
        BadReads.All.Unpair.(chrName).(sampleName{1}) = [];
        BadReads.All.FarPairs.(chrName).(sampleName{1}) = [];
        TotAllLR.(chrName).(sampleName{1}) = [];
        UniqAllLR.(chrName).(sampleName{1}) = [];
        continue
    end
    
    %saving relevant half of reads which didn't map well enough
    BadMapAllData = AllData(AllData(:,3)<=MAPQmin & AllData(:,4)>0,1:4);
    BadReads.All.LowMapQ.(chrName).(sampleName{1}) = [BadMapAllData(:,1) BadMapAllData(:,1) + BadMapAllData(:,4)];
    AllData = AllData(AllData(:,3)>MAPQmin,1:4); %keeping reads which mapped well
    
    % saving unpaired reads and pairs that are too far
    BadReads.All.Unpair.(chrName).(sampleName{1}) = AllData(AllData(:,4)==0,1);
    BadReads.All.FarPairs.(chrName).(sampleName{1}) = [AllData(AllData(:,4)>20000,1) AllData(AllData(:,4)>20000,1) + AllData(AllData(:,4)>20000,4)];
            
    %keep paired reads where left in col 1 and right is not too far
    AllData = AllData(AllData(:,4)>0 & AllData(:,4)<20000,1:4);
    
    % Saving left & right ends (of paired reads)
    TotAllLR.(chrName).(sampleName{1}) = [AllData(:,1) AllData(:,1) + AllData(:,4)];
    
    % Saving unique left & right ends (of paired reads)
    [~,ia,~] = unique(TotAllLR.(chrName).(sampleName{1}),'rows');
    UniqAllData = AllData(ia,:);
    UniqAllLR.(chrName).(sampleName{1}) = [UniqAllData(:,1) UniqAllData(:,1) + UniqAllData(:,4)];
end