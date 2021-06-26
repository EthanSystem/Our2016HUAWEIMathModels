%% init
% clear;
numSamples=1000;
numLocus=9445;
numGene=300;
%% 第一题
%% 数据提取

genotype_origion=importdata('.\OrigionData\genotype.dat');
%% 将数据转换为矩阵

genotype_cell = cell(size(zeros(numSamples,numLocus)));

for i=1:numSamples+1
    str = genotype_origion{i,1};%取出数据中的一列
    aRowOfData = regexp(str,  ' ', 'split');
    [r,c] = size(aRowOfData);
    for j=1:1:c
        genotype_cell{i,j}=aRowOfData{j};
    end
end
%% 建立 genotype 里的位点的序号和 gene_info 里的序号的映射关系

for j=1:1:numLocus
    str010=genotype_cell{1,j};
    str020='rs';
    str030=strrep(str010,str020,'');
    mapOfGenotype_cell{1,j}=num2str(j);
    mapOfGenotype_cell{2,j}=str010;
end
%% 编码成 0 1 2

for j=1:1:9445
    letter(1)='A';letter(2)='T';letter(3)='C';letter(4)='G';
    judge(:,1)=strcmp(genotype_cell(:,j),[letter(1),letter(1)]);
    sum1(1)=sum(judge(:,1));
    judge(:,2)=strcmp(genotype_cell(:,j),[letter(2),letter(2)]);
    sum1(2)=sum(judge(:,2));
    judge(:,3)=strcmp(genotype_cell(:,j),[letter(3),letter(3)]);
    sum1(3)=sum(judge(:,3));
    judge(:,4)=strcmp(genotype_cell(:,j),[letter(4),letter(4)]);
    sum1(4)=sum(judge(:,4));
    [sum1_sorted,sum1_idx]=sort(sum1,'descend');
    goodLetter(1)=letter(sum1_idx(1));
    goodLetter(2)=letter(sum1_idx(2));
    judge5=strcmp(genotype_cell(:,j),[goodLetter(1),goodLetter(2)]);
    sum2=sum(judge5);
    judge6=strcmp(genotype_cell(:,j),[goodLetter(2),goodLetter(1)]);
    sum3=sum(judge6);
    bp{1}=[goodLetter(1),goodLetter(1)];
    bp{3}=[goodLetter(2),goodLetter(2)];
    if sum2>=sum3
        bp{2}=[goodLetter(1),goodLetter(2)];
    else
        bp{2}=[goodLetter(2),goodLetter(1)];
    end
    genotype_cell(:,j)=strrep(genotype_cell(:,j),bp(2),'1');
    genotype_cell(:,j)=strrep(genotype_cell(:,j),bp(1),'0');
    genotype_cell(:,j)=strrep(genotype_cell(:,j),bp(3),'2');
    
end



%% 第二题
% 不在本文件中


%% 第三题
%% 数据提取

for j=1:numGene
    aGenoinfo=eval(['importdata(''.\OrigionData\gene_info\gene_',num2str(j),'.dat'');']);
    for i=1:size(aGenoinfo,1)
        Genes_cell{i,j}=aGenoinfo{i,1};
    end
end
%% 数据映射
%% 生成位点在基因的位置标记


for j=1:numGene
    for i=1:size(Genes_cell,1)
        mark(i,j)=iscellstr(Genes_cell(i,j));
    end
end

%% 生成位点在基因的位置对应的序号值
for j=1:numGene
    for i=1:sum(mark(:,j),1)
        [count,idx]=ismember(mapOfGenotype_cell,Genes_cell{i,j});
        id=find(count(2,:),1);
        indexOfGenotype(i,j)=id;  % index of each genotypes that at each genes in map of genotype_cell.mat.
    end
end

%% 把每个位点所在的基因信息标记在 mapOfGenotype_cell 的第三行
for k=1:1:numLocus
    [row,col]=find(indexOfGenotype==k);
    mapOfGenotype_cell{3,k}=col;
end
clear col;
save('./ResultData/indexOfGenotype.txt','indexOfGenotype');
%% 第三题计算
%% 导入经过MAF算法筛选后的111个位点中的每个位点的序号值数据，
temp=importdata('.\OrigionData\WeiDianIdByMAF.txt');
temp=temp';
numLocus_MAF=size(temp,2);
IndexOfEachGenotype_MAF(1,:)=1:1:numLocus_MAF;
IndexOfEachGenotype_MAF(2,:)=temp(:);


%% 计算每个基因的权值
weightOfEachGene_MAF=zeros(2,numGene);
weightOfEachGene_MAF(1,:)=1:1:numGene;
for j=1:numLocus_MAF
    weightOfEachGene_MAF(2,mapOfGenotype_cell{3,IndexOfEachGenotype_MAF(2,j)})=weightOfEachGene_MAF(2,mapOfGenotype_cell{3,IndexOfEachGenotype_MAF(2,j)})+1;
end
[result , idx]=sort(weightOfEachGene_MAF(2,:),'descend');
numGene_MAF=numGene-sum(weightOfEachGene_MAF_sorted_result==0);
weightOfEachGene_MAF_sorted_result=result(1:numGene_MAF);
weightOfEachGene_MAF_sorted_index=idx(1:numGene_MAF);

%% 绘图
x =1:1:numGene_MAF;
bar(x ,weightOfEachGene_MAF_sorted_result);
xlabel('基因编号（按权值降序排序）');ylabel('基因权值');


%% 附录
%% 计算权值通用方法
weightOfEachGene=zeros(2,numGene);
weightOfEachGene(1,:)=1:1:numGene;

%
for j=1:numLocus_MAF
    weightOfEachGene(2,j)=IndexOfEachGenotype_MAF
end

for j=1:numGene
    for i=1:sum(mark(:,j),1)
        temp=indexOfGenotype(i,j);
        weightOfEachGenotypeAtSameGene=IndexOfEachGenotype_MAF(2,temp);
        weightOfEachGene(2,j)=weightOfEachGene(2,j)+weightOfEachGenotypeAtSameGene;
    end
end
%% 将每个基因的权值降序排序
[weightOfEachGene_sorted_result,weightOfEachGene_sorted_index]=sort(weightOfEachGene(2,:),'descend');

%% 计算贡献率，选取累计贡献率在50之前的所有基因
for j=1:1:numGene
    sumOfSortedWeightInGenes(j)=sum(weightOfEachGene_sorted_result(1:j));
end
for j=1:1:numGene
    contributionRateOfSumOfSortedWeightInGenes(j)=sumOfSortedWeightInGenes(j)/sumOfSortedWeightInGenes(numGene);
end
for j=1:1:numGene
    if contributionRateOfSumOfSortedWeightInGenes(j)>0.5
        index_contributionRateOfSumOfWeightInGenes=j;
        break;
    end
end

contributionGenes=weightOfEachGene_sorted_index(1:index_contributionRateOfSumOfWeightInGenes);


%% 问题4
%% 导入数据
%%
temp041=importdata('.\OrigionData\Q4OrigionData\id\1.txt');
temp042=importdata('.\OrigionData\Q4OrigionData\id\2.txt');
temp043=importdata('.\OrigionData\Q4OrigionData\id\3.txt');
temp044=importdata('.\OrigionData\Q4OrigionData\id\4.txt');
temp045=importdata('.\OrigionData\Q4OrigionData\id\5.txt');
temp046=importdata('.\OrigionData\Q4OrigionData\id\6.txt');
temp047=importdata('.\OrigionData\Q4OrigionData\id\7.txt');
temp048=importdata('.\OrigionData\Q4OrigionData\id\8.txt');
temp049=importdata('.\OrigionData\Q4OrigionData\id\9.txt');
temp0410=importdata('.\OrigionData\Q4OrigionData\id\10.txt');
%%

IndexOfEachGenotype_Q4=zeros(numSamples,10);
%%
IndexOfEachGenotype_Q4(1:size(temp041,1),1)=temp041;
IndexOfEachGenotype_Q4(1:size(temp042,1),2)=temp042;
IndexOfEachGenotype_Q4(1:size(temp043,1),3)=temp043;
IndexOfEachGenotype_Q4(1:size(temp044,1),4)=temp044;
IndexOfEachGenotype_Q4(1:size(temp045,1),5)=temp045;
IndexOfEachGenotype_Q4(1:size(temp046,1),6)=temp046;
IndexOfEachGenotype_Q4(1:size(temp047,1),7)=temp047;
IndexOfEachGenotype_Q4(1:size(temp048,1),8)=temp048;
IndexOfEachGenotype_Q4(1:size(temp049,1),9)=temp049;
IndexOfEachGenotype_Q4(1:size(temp0410,1),10)=temp0410;
% for j=1:1:10
%     eval(['IndexOfEachGenotype_Q4(:,',j,')=temp04',j,';']);
% end
%% 筛选出十个权值里的位点的重复项
mark_Q4=zeros(numSamples,10);
repeatted=zeros(1,numLocus);
for i=1:1:numSamples
    for j=1:1:10
        if IndexOfEachGenotype_Q4(i,j)==0
            break;
        end
        repeatted(1,IndexOfEachGenotype_Q4(i,j))= repeatted(1,IndexOfEachGenotype_Q4(i,j))+1;
    end
end

%% 计算重复项的权值
[repeatted_result,repeatted_index]=sort(repeatted,'descend');

%%
% end