%% create a KEGG construct for the genes that are in p.gene_name(inds)
% NO NEED TO RUN. FILE IS SAVED AND CAN BE LOADED FROM THE DIRECTORY

function [my_kegg,inds_above_th] = produce_kegg_struct_and_plot_violin(gene_names,residual,THRESH,ylab,ytck,ttl)
addpath('X:\Yotam\matlab_projects\spatial_sorting_thesis\1_code\2_functions');
addpath('X:\innaa\MATLAB\');
load('X:\Shalevi\liver_single_cell_rnaseq\bioinformatic_analysis_with_other_databases\Kegg_gene_sets\kegg_pathways.mat');

% L = length(kegg_pathways);
L = 147; % not including the pathlogical keggs

for i = 1:L
    name2 = lower(kegg_pathways(i).name);
    name2(findstr(name2,'_'))=' ';
    name = [upper(name2(1:4)) ': ' upper(name2(6)) lower(name2(7:end))];
    my_kegg(i).name = name;
    my_kegg(i).genes_orig=kegg_pathways(i).genes;
    [b,in] = ismember(lower(kegg_pathways(i).genes),lower(gene_names));
    ind = in(b);
    my_kegg(i).n =      length(ind);
    my_kegg(i).frac =   length(ind)/length(kegg_pathways(i).genes);
    if (~isempty(ind))
        my_kegg(i).ind =         ind;
        my_kegg(i).genes =       gene_names(ind);
        my_kegg(i).resid =       residual(ind);
        my_kegg(i).resid_med =   nanmean(residual(ind));
        my_kegg(i).resid_sem =   std(residual(ind))/sqrt(length(residual(ind)));
        
    end;
end;

% ADD kegg_Carbohydrate_digestion_and_absorption.txt
carbo=read_gene_list('X:\Yotam\gene_sets\kegg_Carbohydrate_digestion_and_absorption.txt');
my_kegg(L+3).name = 'kegg Carbohydrate digestion and absorption';
[b,in]=ismember(lower(carbo),lower(gene_names));
ind=in(b)';
my_kegg(L+3).n=length(ind);
my_kegg(L+3).frac=length(ind)/length(carbo);
if(~isempty(ind))
    my_kegg(L+3).genes=gene_names(ind);
    my_kegg(L+3).ind=ind;
    my_kegg(L+3).resid=residual(ind);
    my_kegg(L+3).resid_med=nanmedian(residual(ind));
end
my_kegg(L+3).genes_orig=carbo;

% plot only pathways which have a reprenation above thresh in the dataset
inds_above_th = [];
for i = 1: length(my_kegg)
    if my_kegg(i).frac >= THRESH & length(my_kegg(i).genes) >= 10 |  ~isempty(find(strcmpi(my_kegg(i).name,'KEGG: Rna polymerase')))
        inds_above_th = [inds_above_th i];
    end
end

keggs2remove = {'KEGG: Oocyte meiosis','KEGG: Cardiac muscle contraction','KEGG: Vascular smooth muscle contraction'};

inds_above_th = setdiff(inds_above_th,find_indices_in_mat({my_kegg.name},keggs2remove));

% uncomment to plot only a set of keggs
% keggs2plot = {'KEGG: Butanoate metabolism','KEGG: Metabolism of xenobiotics by cytochrome p450',...
%     'kEGG: Oxidative phosphorylation','KEGG: Glutathione metabolism','KEGG: Fatty acid metabolism',...
%     'KEGG: Abc transporters','KEGG: Glycerophospholipid metabolism','KEGG: Peroxisome','KEGG: Ribosome',...
%     'KEGG: Spliceosome','KEGG: Citrate cycle tca cycle'}; 
% 
% inds_above_th = find_indices_in_mat({my_kegg.name},keggs2plot);
% %

res=[my_kegg(inds_above_th).resid_med];
[srm,orm]=sort(res,'ascend');

indicator_all = [];
vec_all = [];
for i= 1:length(res)
    vec_all = [vec_all;(my_kegg(inds_above_th(i)).resid)];
    indicator_all=[indicator_all;repmat(i,length(my_kegg(inds_above_th(i)).resid),1)];
end

% produce a violinplot
U=unique(indicator_all);
med=zeros(length(U),1);
for i=1:length(U)
    indd=find(indicator_all==U(i));
    med(i)=median(vec_all(indd));
end
[y,ord]=sort(med,'ascend');
indicator_all2=indicator_all;
for i=1:length(ord)
    indicator_all2(find(indicator_all==ord(i)))=i;
end
addpath('X:\innaa\MATLAB');
figure('units','normalized','outerposition',[0 0 1 1]);
vio = violinplot(vec_all,indicator_all2);
title (ttl);
set(gca,'xtick',[1:1:length(med)] ,'xticklabels',{my_kegg(inds_above_th(ord)).name},'XTickLabelRotation',45);
ylabel(ylab);
grid on;
box on;

cg = [colorGradient([0 .75 0],[1 .75 0],5); colorGradient([1 .75 0],[1 0 0],length(res)-5)];
for i = 1 : length(med), vio(i).ViolinColor = cg(i,:); end
