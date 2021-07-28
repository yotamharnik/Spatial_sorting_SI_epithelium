% load files into worksapce
current_dir = cd;
addpath([current_dir,'\04_matlab_functions\']);

load([current_dir,'\02_processed_data\0_mRNA_clean_UMI_SST_M1-M4.mat']);
load([current_dir,'\02_processed_data\0_villus_zonation_scRNAseq_moor_2018.mat']);
sc.mn_mx = sc.mn./max(sc.mn,[],2);
load([current_dir,'\02_processed_data\1_Protein_clean_iBAQ_SST_M1-M4.mat']);
load([current_dir,'\02_processed_data\2_Protein_mRNA_SC_parsed_SST_M1-M4.mat']);

%% examine mRNA-protein zonation
% calcualte COM
mRNA_com =    calculate_com_mat(sst.sc_mean);
protein_com = calculate_com_mat(sst.protein_norm);

% plot a scatter plot with some integrators and non-integrators proteins(FIGURE 4b)
gg_fit = {'Polr2e','Lypd8','Reg3b','H2-ab1','Cdh1','Ada'};
gg_com = calculate_com_mat(sst.sc_mean(find_indices_in_mat(sst.gene_name,gg_fit),:));
[~,ord] = sort(gg_com);
gg_fit = gg_fit(ord);

gg_integrators = {'Casp6','Nlrp6','Cpt1a','Pck1','Fabp2','slc5a9'};
% gg_com = calculate_com_mat(sst.sc_mn(find_indices_in_mat(sst.gene_name,gg_integrators),:));
% [~,ord] = sort(gg_com);
% gg_integrators = gg_integrators(ord);

[a1,a2]=ismember(lower(gg_fit),lower(sst.gene_name));
ind_genes_fit=a2(a1);

[a1,a2]=ismember(lower(gg_integrators),lower(sst.gene_name));
ind_genes_non_fit=a2(a1);

THRESH = 10^-5;
index = find(max(sst.protein_norm,[],2) > THRESH & max(sst.sc_mean,[],2) > THRESH & sst.protein_mice_count >= 3 & max(sst.protein_cov,[],2) <0.5);

figure;
scatter(mRNA_com(index),protein_com(index),25,[0.55 0.55 0.55],'filled'); hold on;
title('Center of mass (COM) for Protein vs mRNA');
xlabel('mRNA Center of mass');
ylabel('Protein Center of mass');
ylim([0 1]);
set(gca,'xtick',max(xlim),'xticklabel',{'Villus tip'},'ytick',ylim,'yticklabel',{'Villus bottom','Villus tip'});
[r,pp] = corr(mRNA_com(index),protein_com(index),'type','Spearman','rows','complete');
str = ['R_S_p_e_a_r_m_a_n = ',num2str(r,'%.2f')];
text(min(get(gca, 'xlim'))*1.1, max(get(gca, 'ylim'))*0.9,str,'color',[0.6 0.6 0.6]);
plot(xlim, xlim,'k--');
c_fit = [0.9290 0.6940 0.1250];
c_int = [0.6350 0.0780 0.1840];
scatter(mRNA_com(ind_genes_fit),protein_com(ind_genes_fit),75,c_fit,'filled');
text(mRNA_com(ind_genes_fit),protein_com(ind_genes_fit),sst.gene_name(ind_genes_fit),'FontSize',12,'FontWeight','bold');
scatter(mRNA_com(ind_genes_non_fit),protein_com(ind_genes_non_fit),75,c_int,'filled');
text(mRNA_com(ind_genes_non_fit),protein_com(ind_genes_non_fit),sst.gene_name(ind_genes_non_fit),'FontSize',12,'FontWeight','bold');

% axis tight;
box on;
axis square;
grid minor;

% also plot the zonation profiles of these gene (FIGURE 4a 4c)
% gg_integrators = {'Nlrp6','fbp1','fbp2','Pck1','Pck2','slc7a7','slc7a8','Slc15a1','Slc5a1','Slc5a9','slc2a2','cpt1a','Dld','Fabp2','Enpep'};


gg = [gg_integrators,gg_fit];

[m,n] = get_subplot_size(gg);

num_zones = length(zones2include);

figure;
for i=1:length(gg)
    subplot(6,2,i)
    indd = find(strcmpi(sc.gene_name,gg{i}));
    vec    = sc.mn(indd,2:end);
    vec_se = sc.se(indd,2:end);
    plot_patch(linspace(0,1,length(vec)),smoothdata(vec/mean(vec),'loess',7),smoothdata(vec_se/mean(vec),'loess',7),'b');
%     plot_patch(linspace(0,1,length(vec)),vec/max(vec),vec_se/max(vec),'b');
    hold on;
    indd = find(strcmpi(sst.gene_name,gg{i}));
    vec    = sst.protein_norm(indd,:);
    vec_se = sst.protein_sem(indd,:);
    plot_patch(linspace(0,1,length(vec)),smoothdata(vec/mean(vec),'loess',7),smoothdata(vec_se/mean(vec),'loess',7),'r');
%     plot_patch(linspace(0,1,length(vec)),vec/max(vec),vec_se/max(vec),'r');
    hold on;
    set(gca,'xtick',linspace(0,1,size(p(5).protein_norm_mean(:,zones2include),2)));
    if i > 10
        title(gg{i},'Color',c_fit);
    else
        title(gg{i},'Color',c_int);
    end
    xticklabels(p(5).sample_name(1:num_zones))
%     set(gca,'YTick',[0 0.5 1])
    grid on;
    set(gca,'fontsize',10);
    set(gcf,'Position',[680   249   633   729]);
end

%% Examine if there is a shift in the protein COM comapred to mRNA (FIGURE 4d)

THRESH = 10^-5;
index = find(max(sst.sc_mean,[],2)>THRESH & max(sst.protein_norm,[],2)>THRESH & sst.protein_mice_count >=3);

com_rna  = calculate_com_mat(sst.sc_mean(index,:),6);
com_prot = calculate_com_mat(sst.protein_norm(index,:),6);

com_krus = [com_rna ; com_prot];
com_indicator = [repmat(1,length(index),1) ; repmat(2,length(index),1)];
pp=ranksum(com_krus(com_indicator==1),com_krus(com_indicator==2));

binss = 95;
figure;
hp = histogram(com_prot,binss,'FaceColor','r'); hold on;
hp.EdgeAlpha = 0.25;
hp.FaceAlpha = 0.3;
hr = histogram(com_rna,binss,'FaceColor','b');
hr.EdgeAlpha = 0.25;
hr.FaceAlpha = 0.3;
xlabel('Zone')
ylabel('desnity')
title(['histogram of Protein and mRNA COM - #bins ' num2str(binss)]);
% xlim([1 5]);
xline(mean(com_prot),'r--','LineWidth',2);
xline(mean(com_rna),'b--','LineWidth',2);
scatter([prctile(com_prot,25),prctile(com_prot,75)],[min(ylim) min(ylim)],'rx');
scatter([prctile(com_rna,25),prctile(com_rna,75)],[min(ylim) min(ylim)],'bx');
legend('Protein center of mass','mRNA center of mass','Protein mean','mRNA mean','25-75 Protein','25-75 mRNA');
text(min(get(gca, 'xlim'))*1.1, max(get(gca, 'ylim'))*0.9,['p_{ranksum}=' num2str(pp,'%.1e'),' n=' num2str(length(com_krus)/2)]); 
% set(gcf,'Position',[445 554 1156 420]);

% add a density line
y = 0:0.1:6;
mu = mean(com_rna);
sigma = std(com_rna);
f = exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi))*100;
plot(y,f,'LineWidth',2,'Color','b');
hold on;
mu = mean(com_prot);
sigma = std(com_prot);
f = exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi))*100;
plot(y,f,'LineWidth',2,'Color','r');

%% violin plots of mRNA-protein correlations in intestine and liver
THRESH=10^-5;
DYNAMIC_RANGE=1;
corr_all_intestine=NaN*ones(length(sst.gene_name),1);
indin = find(max(sst.protein_norm,[],2) > THRESH );
indin=intersect(indin,find((max(sst.protein_norm,[],2)-min(sst.protein_norm,[],2))./mean(sst.protein_norm,2)>DYNAMIC_RANGE));
for i=1:length(indin),
    corr_all_intestine(indin(i))=corr(sst.protein_norm(indin(i),:)',sst.sc_mean(indin(i),:)','type','spearman');
end

% liver
TT=load('X:\Yotam\Analysis\liver_intestine.mat');
sst_liver=TT.liver;
corr_all_liver=NaN*ones(length(sst_liver.gene_name),1);
indin = find(max(sst_liver.mat_protein,[],2) > THRESH );
indin=intersect(indin,find((max(sst_liver.mat_protein,[],2)-min(sst_liver.mat_protein,[],2))./mean(sst_liver.mat_protein,2)>DYNAMIC_RANGE));
for i=1:length(indin),
    corr_all_liver(indin(i))=corr(sst_liver.mat_protein(indin(i),:)',sst_liver.mat_rna(indin(i),:)','type','spearman');
end

vec=[corr_all_intestine;corr_all_liver];
indicator=[ones(length(corr_all_intestine),1);2*ones(length(corr_all_liver),1);];
figure;
vio = violinplot(vec,indicator);
vio(1).ViolinColor = [0.73 0.73 0];
vio(2).ViolinColor = [0.52 0 0.52];
pp=ranksum(corr_all_liver(~isnan(corr_all_liver)),corr_all_intestine(~isnan(corr_all_intestine)))
ylabel('Correlations between mRNA and proteins');
set(gca,'XTicklabel',{'Intestine','Liver'});
box on;
title(['P=' num2str(pp,'%.1e')]);
line(xlim,[0 0],'color','k','linestyle','--');
set(gca,'fontsize',12);
axis tight

num_neg_liver=100*length(find(corr_all_liver<0))./length(find(~isnan(corr_all_liver)));
num_neg_intestine=100*length(find(corr_all_intestine<0))./length(find(~isnan(corr_all_intestine)));
text(1.8,-0.2, [num2str(num_neg_liver,'%.1f') '%']);
text(0.8,-0.2, [num2str(num_neg_intestine,'%.1f') '%']);
ylim([-1.1 1.1]);