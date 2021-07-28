% add the path to all the necessary functions and defince working directory
current_dir = cd;
addpath([current_dir,'\04_matlab_functions\']);

%% process raw scRNAseq data from 'Moor et al cell 2018' [max Normalzied UMI count]
% LOAD AND ORGENIZE RAW DATA

moor_table = readtable([current_dir,'\01_RNAseq_msProteomics_flies\Villus_zonation_Cell_version_August_2018.xlsx']);

sc.gene_name = moor_table.GeneName;
sc.mn = [moor_table.Crypt_mean moor_table.V1_mean moor_table.V2_mean moor_table.V3_mean moor_table.V4_mean moor_table.V5_mean moor_table.V6_mean];
sc.se = [moor_table.Crypt_SE   moor_table.V1_SE   moor_table.V2_SE   moor_table.V3_SE   moor_table.V4_SE   moor_table.V5_SE   moor_table.V6_SE];

%% save the orgenzied scRNAseq data

folder_name = '\02_processed_data';
file_name = '\0_villus_zonation_scRNAseq_moor_2018.mat';

% save([current_dir,folder_name,file_name],'sc');

%% process raw mRNA data [UMI count]
% LOAD AND ORGENIZE RAW DATA

mn = readtable([current_dir,'\01_RNAseq_msProteomics_flies\zUMIs_Table_SST_WT_M1-M4.txt']);

st.gene_name    = table2array(mn(:,2));
st.sample_name  = mn.Properties.VariableNames(3:end);
st.zone_name    = {'V1','V2','V3','V4','V5','V6'};
st.mice_name    = {'M1','M2','M3','M4'};
st.mat          = table2array(mn(:,3:end));

% define zone and mice flags
zones2include = [1 2 3 4 5 6];
num_zones = length(zones2include);

mice2include    = [1 2 3 4];
num_mice  = length(mice2include);


%% CLEAN RAW DATA

% remove abosrptive markers
% find markers from Monaco et. al 2020
goblet_genes_inds  = find_indices_in_mat(st.gene_name,read_gene_list([current_dir,'\01_RNAseq_msProteomics_flies\markers_goblet.txt']));
EEC_genes_inds     = find_indices_in_mat(st.gene_name,read_gene_list([current_dir,'\01_RNAseq_msProteomics_flies\markers_eec.txt']));
tuft_genes_inds    = find_indices_in_mat(st.gene_name,read_gene_list([current_dir,'\01_RNAseq_msProteomics_flies\markers_tuft.txt']));
Paneth_genes_inds = find_indices_in_mat(st.gene_name,read_gene_list([current_dir,'\01_RNAseq_msProteomics_flies\markers_paneth.txt']));
indout          = unique([goblet_genes_inds,EEC_genes_inds,tuft_genes_inds,Paneth_genes_inds]);
% find duplicates
[~,ind_unique,~] = unique(st.gene_name);

% reomve duplicate genes and markers from data
indin           = intersect(setdiff(1:length(st.gene_name),indout),ind_unique);

st.gene_name    = st.gene_name(indin);
st.mat          = st.mat(indin,:);
st.m1           = st.mat(:,1:6);
st.m2           = st.mat(:,7:12);
st.m3           = st.mat(:,13:18);
st.m4           = st.mat(:,19:24);

%% NORMALISE AND CALCULATE MEAN AND S.E.M

% normalize to the sum
st.mat_norm     = st.mat./nansum(st.mat);
st.m1_norm      = st.mat_norm(:,1:6);
st.m2_norm      = st.mat_norm(:,7:12);
st.m3_norm      = st.mat_norm(:,13:18);
st.m4_norm      = st.mat_norm(:,19:24);

% alucate space to hold the normlized means across mice
st.mean         = zeros(size(st.mat,1),num_zones);
st.sem          = zeros(size(st.mat,1),num_zones);
st.mean_norm    = zeros(size(st.mat,1),num_zones);
st.sem_norm     = zeros(size(st.mat,1),num_zones);


for i=1:size(st.mat_norm,1)
    % extract the data
    vec(1,:)      = st.m1(i,zones2include);
    vec(2,:)      = st.m2(i,zones2include);
    vec(3,:)      = st.m3(i,zones2include);
    vec(4,:)      = st.m4(i,zones2include);
    
    vec_norm(1,:) = st.m1_norm(i,zones2include);
    vec_norm(2,:) = st.m2_norm(i,zones2include);
    vec_norm(3,:) = st.m3_norm(i,zones2include);
    vec_norm(4,:) = st.m4_norm(i,zones2include);
    
    % calcualte mean
    st.mean(i,:)  = nanmean(vec);
    st.sem(i,:)   = nanstd(vec)/sqrt(length(mice2include));
    st.mean_norm(i,:) = nanmean(vec_norm);
    st.sem_norm(i,:)  = nanstd(vec_norm)/sqrt(length(mice2include));
    clear vec
    clear vec_norm
end

%% PREFORM QC

% examine the sum of UMIs per samples
figure;bar(sum(st.mat));
title('Sum of zUMIs per sample','Fontsize',14);
ylabel('Sum of zUMIs','Fontsize',12);
xlabel('sample name','Fontsize',12);
xticks(1:length(st.sample_name));
xticklabels(st.sample_name);
xtickangle(45);

% check if differnet zones have different mRNA extraction rates
sum_all=sum(st.mat);
indicator=[1:num_zones 1:num_zones 1:num_zones 1:num_zones];
kruskalwallis(sum_all,indicator);
xticks(linspace(1,num_zones,num_zones));
xticklabels(st.zone_name);

% check if differnet mice had different mRNA extraction rates
sum_all=sum(st.mat);
indicator=[repmat(1,1,num_zones) repmat(2,1,num_zones)...
           repmat(3,1,num_zones) repmat(4,1,num_zones)];
kruskalwallis(sum_all([1:num_zones,(1:num_zones)+num_zones,...
    (1:num_zones)+num_zones*2,(1:num_zones)+num_zones*3]),indicator);
xticks(1:num_mice);
xticklabels(st.mice_name);

%% save the clean UMI talbe

folder_name = '\02_processed_data';
file_name = '\0_mRNA_clean_UMI_SST_M1-M4.mat';

% save([current_dir,folder_name,file_name],'st');

%% process raw Protein data [proteinGroups intestines]:
% LOAD AND ORGENIZE RAW DATA 

% load the data from the xls file to a struct named p
path_to_flies = [current_dir,'\01_RNAseq_msProteomics_flies\'];

% grab all xlsx flies in the raw data folder and orgenize each mice in the p struct
listing = dir(path_to_flies);
listing = {listing(find(contains({listing.name},'SI_11885'))).name};

for i = 1 : length(listing)
    [A,B] = xlsread([path_to_flies listing{i}]);
    p(i).ID =          B(2:end,1);
    p(i).sample_name = {['M' num2str(i) '_V1'],['M' num2str(i) '_V2'],['M' num2str(i) '_V3'],...
                        ['M' num2str(i) '_V4'],['M' num2str(i) '_V5'],['M' num2str(i) '_V6']};
    p(i).gene_name =   B(2:end,4);
    p(i).peptide =     A(:,1);
    p(i).raw_iBAQ =    A(:,9:14);
    p(i).mol_weight =  A(:,16);
    clear A B;
end

%% CLEAN RAW DATA

% remove entries with all missing values
for j = 1:length(listing)
    indout = find(sum(p(j).raw_iBAQ,2) == 0);
    indout_missing{j} = indout;
end

% find most highly exprsse protein at the RNA level for protein with several proteins in the same gene_name
load('X:\Yotam\matlab_projects\spatial_sorting_thesis\0_data\1_tidy_data\0_mRNA_processed_data_SST_WT_M1-M4.mat');
p_orig =  p;
for j = 1:length(listing)    
    multi_proteins_inds = find(contains(p(j).gene_name,';'));
    for i = 1 : length(multi_proteins_inds)
        [proteins] = split(p(j).gene_name(multi_proteins_inds(i)),';');
        high_protein_ind = find_most_high_exp_protein_in_mrna(proteins,st);
        p(j).gene_name(multi_proteins_inds(i)) = proteins(high_protein_ind);
    end
end

% find lab contaminents based on CON_ or missing gene names
for j = 1:length(listing)
    indout = [];
    for i = 1 : length(p(j).ID)
        if p(j).gene_name(i) == ""
            indout=[indout i];
        end
    end
    indout_con{j} = indout;
end

% find proteins who have zreo expression in the single-cell gene expression data
for j = 1:length(listing)
    indout = [];
    for i = 1 : length(p(j).gene_name)
        index = find(strcmpi(sc.gene_name,(p(j).gene_name(i))));
        sc_sum = sum(sc.mn(index,:));
        if isempty(index) | sc_sum == 0 
            indout=[indout i];
        end
    end    
    indout_con_cs{j} = indout;
end

clear goblet_genes_inds EEC_genes_inds tuft_genes_inds Paneth_genes_inds
% find abosprtive markers from Monaco et al 2020
for j = 1 : length(listing)
    goblet_genes_inds{j}  = find_indices_in_mat(p(j).gene_name,read_gene_list([current_dir,'\01_RNAseq_msProteomics_flies\markers_goblet.txt']));
    EEC_genes_inds{j}     = find_indices_in_mat(p(j).gene_name,read_gene_list([current_dir,'\01_RNAseq_msProteomics_flies\markers_eec.txt']));
    tuft_genes_inds{j}    = find_indices_in_mat(p(j).gene_name,read_gene_list([current_dir,'\01_RNAseq_msProteomics_flies\markers_tuft.txt']));
    Paneth_genes_inds{j}  = find_indices_in_mat(p(j).gene_name,read_gene_list([current_dir,'\01_RNAseq_msProteomics_flies\markers_paneth.txt']));
end

% remove contaminates 
for j = 1 : length(listing)
    clear indout
    indout =        [indout_missing{j}' ,goblet_genes_inds{j}, EEC_genes_inds{j}, tuft_genes_inds{j},Paneth_genes_inds{j}, indout_con{j}, indout_con_cs{j}];
    indin =         setdiff(1:length(p(j).gene_name),unique(indout));
    p(j).gene_name =   p(j).gene_name(indin);
    p(j).peptide =     p(j).peptide(indin);
    p(j).mol_weight =  p(j).mol_weight(indin);
    p(j).raw_iBAQ =    p(j).raw_iBAQ(indin,:);
end


%% NORMALISE

% Assign all mice to the common-proteins list (add missing values)
p_orig = p;

all_genes = [];
for j = 1:num_mice
    all_genes = union(all_genes,p(j).gene_name);
end

for j = 1 : length(listing)    
p(j).gene_name =     all_genes;
p(j).ID =            {};
p(j).raw_iBAQ =      zeros(length(all_genes),length(p(j).sample_name));
p(j).peptide =       zeros(length(all_genes),1);
p(j).mol_weight =    zeros(length(all_genes),1);
end

for i=1:length(all_genes)
    ind1=find(strcmpi(p_orig(1).gene_name,all_genes{i}));
    ind2=find(strcmpi(p_orig(2).gene_name,all_genes{i}));
    ind3=find(strcmpi(p_orig(3).gene_name,all_genes{i}));
    ind4=find(strcmpi(p_orig(4).gene_name,all_genes{i}));
    if ~isempty(ind1)
        p(1).ID(i) =            p_orig(1).ID(ind1(1));
        p(1).raw_iBAQ(i,:) =    p_orig(1).raw_iBAQ(ind1(1),:);
        p(1).peptide(i) =       p_orig(1).peptide(ind1(1));
        p(1).mol_weight(i) =    p_orig(1).mol_weight(ind1(1));
    end
    if ~isempty(ind2)
        p(2).ID(i) =            p_orig(2).ID(ind2(1));        
        p(2).raw_iBAQ(i,:) =    p_orig(2).raw_iBAQ(ind2(1),:);
        p(2).peptide(i) =       p_orig(2).peptide(ind2(1));
        p(2).mol_weight(i) =    p_orig(2).mol_weight(ind2(1));
    end
    if ~isempty(ind3)
        p(3).ID(i) =            p_orig(3).ID(ind3(1));
        p(3).raw_iBAQ(i,:) =    p_orig(3).raw_iBAQ(ind3(1),:);
        p(3).peptide(i) =       p_orig(3).peptide(ind3(1));
        p(3).mol_weight(i) =    p_orig(3).mol_weight(ind3(1));
    end
    if ~isempty(ind4)
        p(4).ID(i) =            p_orig(4).ID(ind4(1));
        p(4).raw_iBAQ(i,:) =    p_orig(4).raw_iBAQ(ind4(1),:);
        p(4).peptide(i) =       p_orig(4).peptide(ind4(1));
        p(4).mol_weight(i) =    p_orig(4).mol_weight(ind4(1));
    end
end

% normalize by sum
for j = 1:length(listing)
    p(j).iBAQ_norm =   p(j).raw_iBAQ./sum(p(j).raw_iBAQ);
%     p(j).iBAQ_norm(find(isnan(p(j).iBAQ_norm)))  = 0 ; 
end

%% PREFORM QC

% examine the sum of iBAQ per samples
figure;
bar([sum(p(1).raw_iBAQ),sum(p(2).raw_iBAQ),sum(p(3).raw_iBAQ),sum(p(4).raw_iBAQ)]);
title('Sum of intensitis per sample','Fontsize',14);
ylabel('Sum of intensitis','Fontsize',12);
xlabel('sample name','Fontsize',12);
xticks(1:num_mice*num_zones);
xticklabels([p(1).sample_name p(2).sample_name p(3).sample_name p(4).sample_name]);
xtickangle(45);


% check if differnet zones have different protein extracion rates
sum_all=[sum(p(1).raw_iBAQ),sum(p(2).raw_iBAQ),sum(p(3).raw_iBAQ),sum(p(4).raw_iBAQ)];
indicator=[1:num_zones 1:num_zones 1:num_zones 1:num_zones];
kruskalwallis(sum_all,indicator);
xticks(linspace(1,num_zones,num_zones));
title('sum of protein intestines per zone')

% check if differnet mice had different protein extraction rates
sum_all=[sum(p(1).raw_iBAQ),sum(p(2).raw_iBAQ),sum(p(3).raw_iBAQ),sum(p(4).raw_iBAQ)];
indicator=[repmat(1,1,num_zones) repmat(2,1,num_zones) repmat(3,1,num_zones) repmat(4,1,num_zones)];
kruskalwallis(sum_all([1:num_zones,(1:num_zones)+num_zones,...
    (1:num_zones)+num_zones*2,(1:num_zones)+num_zones*3]),indicator);
xticks(1:4);
xticklabels(st.mice_name);

%% CALCULATE MEAN AND S.E.M

% make a list of all of the proteins measuerd in the 5 mice.
p(5).gene_name =              all_genes;
p(5).sample_name =            {'V1','V2','V3','V4','V5','V6'};
p(5).protein_norm_mean =      zeros(length(p(5).gene_name),length(p(5).sample_name));
p(5).protein_norm_se =        zeros(length(p(5).gene_name),length(p(5).sample_name));
p(5).protein_norm_median =    zeros(length(p(5).gene_name),length(p(5).sample_name));
p(5).num_of_expressing_mice = zeros(length(p(5).gene_name),1);
p(5).mean_raw_iBAQ =          zeros(length(p(5).gene_name),length(p(5).sample_name));  

for i=1:length(p(5).gene_name)
    clear tab_protein;
    clear tab_protein_raw;
    clear mice_count;
    
    tab_protein = [p(1).iBAQ_norm(i,:); p(2).iBAQ_norm(i,:); p(3).iBAQ_norm(i,:); p(4).iBAQ_norm(i,:)]; 
    tab_protein_raw = [p(1).raw_iBAQ(i,:); p(2).raw_iBAQ(i,:); p(3).raw_iBAQ(i,:); p(4).raw_iBAQ(i,:)]; 
    mice_count = sum([max(p(1).iBAQ_norm(i,:))>0 max(p(2).iBAQ_norm(i,:))>0 max(p(3).iBAQ_norm(i,:))>0 max(p(4).iBAQ_norm(i,:))>0]);

    p(5).protein_norm_mean(i,:) =      mean(tab_protein,1);
    p(5).protein_norm_median(i,:) =    median(tab_protein,1);
    p(5).protein_norm_se(i,:) =        std(tab_protein,[],1)/sqrt(size(tab_protein,1));
    p(5).num_of_expressing_mice(i) =   mice_count;
    p(5).mean_raw_iBAQ(i,:) =          mean(tab_protein_raw,1);
end



%% save the clean iBAQ talbe

folder_name = '\02_processed_data';
file_name = '\1_Protein_clean_iBAQ_SST_M1-M4.mat';

% save([current_dir,folder_name,file_name],'p');

%% integrate data sets: Spatial sorting bulk-mRNAseq, Spatial sorting mass-spec Proteomics and singe-cell RNAseq spatial reconstruction
cross_genes = intersect(lower(p(5).gene_name),lower(sc.gene_name));

sst.gene_name           = cross_genes;

sst.protein_norm        = zeros(length(cross_genes),num_zones);
sst.protein_sem         = zeros(length(cross_genes),num_zones);
sst.protein_cov         = zeros(length(cross_genes),num_zones);
sst.protein_mice_count  = zeros(length(cross_genes),1);

sst.mRNA_norm           = zeros(length(cross_genes),num_zones);
sst.mRNA_sem            = zeros(length(cross_genes),num_zones);

sst.raw_iBAQ            = zeros(length(cross_genes),num_zones);
sst.raw_umiz            = zeros(length(cross_genes),num_zones);

sst.sc_mean             = zeros(length(cross_genes),num_zones);
sst.sc_sem              = zeros(length(cross_genes),num_zones);

for i = 1 : length(sst.gene_name)
    indp  = find(strcmpi( p(5).gene_name,  sst.gene_name{i}));
    indr  = find(strcmpi( st.gene_name,    sst.gene_name{i}));
    indsc = find(strcmpi( sc.gene_name,    sst.gene_name{i}));
    if ~isempty(indp)
        indp = indp(1);
        sst.protein_norm(i,:)           = p(5).protein_norm_mean(indp,zones2include);
        sst.protein_norm_median(i,:)    = p(5).protein_norm_median(indp,zones2include);
        sst.protein_sem(i,:)            = p(5).protein_norm_se(indp,zones2include);
        sst.protein_cov(i,:)            = p(5).protein_norm_se(indp,zones2include)./p(5).protein_norm_mean(indp,zones2include);        
        sst.raw_iBAQ(i,:)               = p(5).mean_raw_iBAQ(indp,zones2include);
        sst.protein_mice_count(i)       = p(5).num_of_expressing_mice(indp);
    end
    if ~isempty(indr)
        indr = indr(1);
        sst.mRNA_norm(i,:)              = st.mean_norm(indr,zones2include);
        sst.mRNA_sem(i,:)               = st.sem_norm(indr,zones2include);
        sst.raw_umiz(i,:)               = st.mean(indr,zones2include);
    end
    if ~isempty(indsc)
        sst.sc_mean(i,:)                = sc.mn(indsc,2:end);
        sst.sc_sem(i,:)                 = sc.se(indsc,2:end);
    end
end

% save the parsed sst structer
folder_name = '\02_processed_data';
file_name = '\2_Protein_mRNA_SC_parsed_SST_M1-M4.mat';

% save([current_dir,folder_name,file_name],'sst');
