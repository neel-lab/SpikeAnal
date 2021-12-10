function mutation=spike_tall (varargin)
% This MATLAB script is used to read Spike mutation fasta files from GIASID and
% quantify mutations at each site, with focus on N- and O-linked
% glycosylation site mutations. But, of couse, it can be used for other
% site mutation quantitation also

% Sriram Neelamegham 11/16/2021

tic
clear;
if nargin==1
    inputFile=varargin{1};
else
    inputFile='spikeprot1113.txt';
end
ref='MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQGVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT*';
ds = tabularTextDatastore(inputFile,'NumHeaderLines',0);  % input file from GISAID, only changed file name to .txt rather than .fasta due to MATLAB quirk
T=tall(ds);             % read data into tall array
szT=gather(size(T,1));  % get size of tall array
n=500000;               % analyze 500,000 reads at a time
len=1270;               % include up to 4 deletions (1274-1270)
nX=3;                   % include up to 3X in read
nblank=3;               % maximum of 3 insertions
idx=floor(szT/n);       % number of cycles for the program
totalReads=0;
mutation=zeros(1,length(ref));  % initialize final result array

for i=1:idx
    Spike=table2cell(gather(T((i-1)*n+1:i*n,:)));  % read data
    Spike(contains(Spike,'>'))=[];                 % remove fasta header
    Spike(cellfun(@(x) length(x)<len, Spike))=[];  % remove short stuff
    Spike(cell2mat(cellfun(@(x) length(find(x=='X')), Spike, 'UniformOutput', false))>nX)=[]; % remove if >nX
    [uSpike,~,ic]=unique(Spike);                  % find unique seq. and store as uSpike
    uCt=histc(ic,unique(ic));                      % calculate freq. of unique seq. in uSpike
    totalReads=totalReads+sum(uCt);
    [~,AuSpike]=cellfun(@(x) nwalign(x,ref), uSpike, 'UniformOutput', false); % align uSpike against ref.
    
    % this loop replaces all X with the most likely amino acid based on reference sequence
    ans1=cellfun(@(x) strfind(x(1,:),'X'), AuSpike,'UniformOutput', false);   % mark position with 'X' in each AuSpike first row
    tf=~cellfun(@(x) isempty(x), ans1);    % which AuSpike have X
    for j=1:size(ans1,1)
        if tf(j)
            idx1=cell2mat(ans1(j));
            idx1Rep=(idx1-1)*3+1;
            for m=1:length(idx1Rep)
                AuSpike{j}(idx1Rep)=AuSpike{j}(idx1Rep+2);
                AuSpike{j}(idx1Rep+1)='|';
            end
        end
    end
    
    % this loop is to correct sequence position for insertion mutations
    ans1=cellfun(@(x) strfind(x(3,:),'-'), AuSpike,'UniformOutput', false);   % mark position with '-' in reference as this corresponds to insertions
    tf=~cellfun(@(x) isempty(x), ans1);    % which AuSpike have '-'
    for j=size(ans1,1):-1:1
        if tf(j)
            idx1=cell2mat(ans1(j));
            idx1Rep=(idx1-1)*3+1;
            m=length(idx1Rep);
            while (m>0)
                AuSpike{j}(idx1Rep:idx1Rep+2)=[];
                AuSpike{j}(idx1Rep)=':';
                sz=length(AuSpike{j})/3;
                if (floor(sz)-sz==0)
                    AuSpike{j,1}=reshape(AuSpike{j,1},[3,sz]);
                else
                    AuSpike(j)=[];
                    uCt(j)=[];
                    m=1;
                end
                m=m-1;
            end
        end
    end
    
    PAuSpike=cellfun(@(x) x(2,:), AuSpike, 'UniformOutput', false); % pick just the Aligned Spike seq. and analyze them for mutation freq.
    parfor k=1:length(ref)
        an=cellfun(@(v)v(k),PAuSpike);
        temp(k)=sum((an~='|').*uCt);
    end
    mutation=mutation+temp;
end
normMutation=mutation/totalReads*100;  % normalized mutation rates based on number of sequences analyzed

% focus on N-glycosylation
idx=[17,61,74,122,149,165,234,282,331,343,603,616,657,709,717,801,1074,1098,1134,1158,1173,1194];
normMutation(idx)
normMutation(idx+3)
totalReads  % total number of reads analayzed
toc
end