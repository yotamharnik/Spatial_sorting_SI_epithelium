close all
clear all
%% Compare spatial covrage of Spatial sorting and scRNAseq reconstruction 
current_dir = cd;
addpath([current_dir,'\04_matlab_functions\']);

%% examine the RNA zonation corrolation between scRNAseq data and spatial sorting

% load the two data sets
load([current_dir,'\02_processed_data\0_mRNA_clean_UMI_SST_M1-M4.mat']);
load([current_dir,'\02_processed_data\0_villus_zonation_scRNAseq_moor_2018.mat']);
sc.mn_mx = sc.mn./max(sc.mn,[],2);

% define which gates to use
st_zones2include = 1:6;
sc_zones2include = 1:7;

% calcualte center of mass (COM) for SC and SST (takes some time)
st_sc_cross_genes.genes    = intersect(lower(sc.gene_name),lower(st.gene_name));
st_sc_cross_genes.st_ind   = zeros(length(st_sc_cross_genes.genes),1);
st_sc_cross_genes.sc_ind   = zeros(length(st_sc_cross_genes.genes),1);
st_com  =  zeros(length(st_sc_cross_genes.genes),1);
sc_com  =  zeros(length(st_sc_cross_genes.genes),1);

for i=1:length(st_sc_cross_genes.genes)
    stind = find(strcmpi(st.gene_name,st_sc_cross_genes.genes{i}));
    st_sc_cross_genes.st_ind(i) = stind(1);
    st_com(i) = nansum(st_zones2include.*st.mean_norm(stind(1),st_zones2include))/nansum(st.mean_norm(stind(1),st_zones2include));
    
    scind = find(strcmpi(sc.gene_name,st_sc_cross_genes.genes{i}));
    st_sc_cross_genes.sc_ind(i) = scind(1);
    sc_com(i)  = nansum(sc_zones2include.*sc.mn(scind(1),sc_zones2include))/nansum(sc.mn(scind(1),sc_zones2include));
end

% plot the corrolation (FIGURE 3F)

THRESH = 10^-5;
inds   = find(max(st.mean_norm(st_sc_cross_genes.st_ind),[],2) > THRESH & max(sc.mn(st_sc_cross_genes.sc_ind),[],2) > THRESH ...
    & ~isnan(st_com) & ~isnan(sc_com));

figure;
xx = st_com(inds);
sort_ind = sc_com(inds);
scatter(xx,sort_ind,10,'filled'); 
xlabel('Spatial sorting bulk RNAseq [Center of mass]');
ylabel('single-cell RNAseq reconstruction [Center of mass]');
xticks([min(xlim) max(xlim)])
xticklabels({'Villus bottom','Villus tip'});
yticks([min(ylim) max(ylim)])
yticklabels({'','Villus tip'});
[r,zonation_pval] = corr(xx,sort_ind,'type','Spearman');
hold on;
text(0.05,0.95,['R_{Spearman}=',num2str(r,'%.2f')],'units','normalized','FontSize',12);
[pp,S] = polyfit(xx,sort_ind,1);
plot([min(xx) max(xx)],[min(pp(1)*xx+pp(2)),max(pp(1)*xx+pp(2))],'k--','LineWidth',1.2);
axis tight;
box on;
axis square;

genes2plot = {'Rpl3','CD9','Nt5c','H2-ab1','Slc15a1','Hist1h1c','Fos','Slc28a2','Nt5e'};

gg_inds = find_indices_in_mat(st_sc_cross_genes.genes,genes2plot);
gg_inds = intersect(gg_inds,inds);
scatter(st_com(gg_inds),sc_com(gg_inds),'ko','filled');

text(st_com(gg_inds),sc_com(gg_inds),st_sc_cross_genes.genes(gg_inds),'FontSize',12);

xticks([2.2 4.6]);
yticks([1 6.6]);
xticklabels({'Villus Bottom','Villus tip'});
yticklabels({'Villus Bottom','Villus tip'});
ytickangle(90);

% plot  profles of scRNAseq data and spatial sorting (FIGURE 3G)
gg_inds = find_indices_in_mat(sc.gene_name,genes2plot);

figure;
for i = 1:length(genes2plot)
    subplot(3,3,i)
    ind_st    = find(strcmpi(st.gene_name,genes2plot{i}));
    ind_sc = find(strcmpi(sc.gene_name,genes2plot{i}));
    if ~isempty(ind_st) & ~isempty(ind_sc)
        vec = st.mean_norm(ind_st,st_zones2include);
        vec_se = st.sem_norm(ind_st,st_zones2include);
        plot_patch(linspace(0,1,length(vec)),vec/max(vec),vec_se/max(vec),'b'); hold on;
        vec = sc.mn(ind_sc,sc_zones2include);
        vec_se = sc.se(ind_sc,sc_zones2include);
        plot_patch(linspace(0,1,length(vec)),vec/max(vec),vec_se/max(vec),[1 .75 0]);
    end 
    ylim([0 max(ylim)]);
    xticks(linspace(0,1,length(st_zones2include)));
    xticklabels(st.zone_name(1:length(st_zones2include)));
    yticks(linspace(0,2,5));
    set(gca,'FontSize',10);
    title(genes2plot{i},'FontSize',20);
    box on;
    grid on;
    if i ==1
        legend('Spatial sorting','','Single-cell reconstruction','Orientation','horizontal');
    end
end
set(gcf,'Position',[730   447   454   520]);

%% Calculate protein zonation

% load the protein data
load([current_dir,'\02_processed_data\1_Protein_clean_iBAQ_SST_M1-M4.mat']);
zones2include = 1:6;

% calcualte Center of mass
p(5).com =    zeros(length(p(5).protein_norm_mean),1);
for i=1:length(p(5).com)
    % scale protein (only divide by max)
    vec = p(5).protein_norm_mean(i,zones2include)./(max(p(5).protein_norm_mean(i,zones2include)));
    p(5).com(i) = sum(linspace(0,1,length(vec)).*vec)/sum(vec);
end

THRESH = 10^-6;
indin = find(max(p(5).protein_norm_mean,[],2)>THRESH & p(5).num_of_expressing_mice >= 2) ;
[y,ord] = sort(p(5).com(indin));

% evaluate zonation statisticly - extract pval for each protein over all zones
zonation_pval = zeros(size(p(5).gene_name(indin),1),1);
indicator = [1:6,1:6,1:6,1:6];

for i=1:size(p(5).gene_name(indin),1)
    if sum(p(5).protein_norm_mean(indin(i),zones2include))~=0
        krus_vec = [p(1).iBAQ_norm(indin(i),zones2include)/mean(p(1).iBAQ_norm(indin(i),zones2include),2) ...
            p(2).iBAQ_norm(indin(i),zones2include)/mean(p(2).iBAQ_norm(indin(i),zones2include),2) ...
            p(3).iBAQ_norm(indin(i),zones2include)/mean(p(3).iBAQ_norm(indin(i),zones2include),2) ...
            p(4).iBAQ_norm(indin(i),zones2include)/mean(p(4).iBAQ_norm(indin(i),zones2include),2)];
        zonation_pval(i) = anovan(krus_vec,{indicator},'display','off');
    end
end

% compute qval
qval = NaN*ones(length(zonation_pval),1);
qval = mafdr(zonation_pval,'BHFDR','True');

% filter genes above qval thresh and above above exp thresh
QVAL_THRESH=0.25;
indin2=intersect(indin,find(qval<QVAL_THRESH));

% plot heatmap of protein zonation (Figure 4A)
genes2plot =  {'Polr2c','SRPR','hspa5','Lypd8','H2-ab1','Pigr','Mdh2','Ndufa5','Acbd5','cyp4b1','Apoa4','Ada'};

gg_com = calculate_com_mat(p(5).protein_norm_mean(find_indices_in_mat(p(5).gene_name,genes2plot),:));
[~,ord_gg] = sort(gg_com);
genes2plot = genes2plot(ord_gg);

ind_genes = find_indices_in_mat(p(5).gene_name(indin(ord)),genes2plot);
[sort_ind,order_genes2plot] = sort(ind_genes);

figure;
protein_carpet = p(5).protein_norm_mean(indin(ord),zones2include)./repmat(max(p(5).protein_norm_mean(indin(ord),zones2include),[],2),1,size(p(5).protein_norm_mean(indin(ord),zones2include),2));
im = imagesc(smoothdata(protein_carpet,'lowess',5));
for i=0:size(protein_carpet,2), line([i+0.5 i+0.5],ylim,'color','k','linewidth',1);end
set(gca,'fontsize',15);
title(['n = ',num2str(length(indin)),' proteins. expression thresh = ',num2str(THRESH),newline,...
    num2str(length(indin2)/length(indin)*100,'%.0f'),'% proteins are significantly zonated qval < ',num2str(QVAL_THRESH)]);
set(gca,'YTick',sort_ind,'YTickLabel',p(5).gene_name(indin(ord(ind_genes(order_genes2plot)))),'FontSize',8);
set(gca,'xtick',1:size(p(5).protein_norm_mean(:,zones2include),2));
xlabel('Facs Gates','fontsize',16);
ylabel('Protein Expression (Normalized by the Maximum)','fontsize',16);
set(gcf,'position',[440    42   411   954]);
c = colorbar ; c.Label.String = 'Relative Epression'; c.FontSize = 10; box on;
set(gca,'xticklabel',{'V1','V2','V3','V4','V5','V6'});

%% plot protein zonation profiles

genes2plot =  {'Polr2c','SRPR','hspa5','Lypd8','H2-ab1','Pigr','Mdh2','Ndufa5','Acbd5','cyp4b1','Apoa4','Ada'};
figure;
for i = 1:length(genes2plot)
    subplot(6,2,i)
    ind    = find(strcmpi(p(5).gene_name,genes2plot{i}));
    if ~isempty(ind) 
        vec = p(5).protein_norm_mean(ind,zones2include);
        vec_se =  p(5).protein_norm_se(ind,zones2include);
%         plot_patch(linspace(0,1,length(vec)),smoothdata(vec,'lowess',5),smoothdata(vec_se,'lowess',5),'k'); hold on;
        plot_patch(linspace(0,1,length(vec)),vec,vec_se,'k'); hold on;

    end
%     ylim([0 max(ylim)]);
    xticks(linspace(0,1,length(zones2include)));
    xticklabels(st.zone_name(1:length(zones2include)));
%     yticks(linspace(0,2,5));
    set(gca,'FontSize',10);
    title(genes2plot{i},'FontSize',12);
    box on;
    grid on;
    ylim([0 max(ylim)]);
end
set(gcf,'Position',[730    42   438   954]);

%% calcualte KEGG pathways zonation (FIGURE 3b)

% load the integrated data
load([current_dir,'\02_processed_data\2_Protein_mRNA_SC_parsed_SST_M1-M4.mat']);

protein_com = calculate_com_mat(sst.protein_norm);

THRESH = 10^-6;
indin = find(max(sst.protein_norm,[],2) > THRESH & sst.protein_mice_count >= 2);

% take only zonmated proteins
sst.protein_max_dynamic_range = (max(sst.protein_norm,[],2) - min(sst.protein_norm,[],2))./mean(sst.protein_norm,2);
THRESH_dynamic = 0.5;
indin = intersect(indin,find(sst.protein_max_dynamic_range >= THRESH_dynamic));

FRAC_THRESH = 0.2;
ylab = 'pathway mean Protein COM' ;
ytck = {'villus bottom','villus tip'};
ttl = '';

[KEGG_pathways,inds_above_th] = produce_kegg_struct_and_plot_violin(sst.gene_name(indin),protein_com(indin),FRAC_THRESH,ylab,ytck,ttl);
ylim([0.2 0.8]);

