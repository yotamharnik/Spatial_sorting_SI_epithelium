function [yyr_mat,yyp_mat] = plot_protein-mrna_zonation_with_delta(gg,NORM,plot_flag);

if nargin <2
    NORM = 0;
    plot_flag = 1;
end
if nargin <3
    plot_flag = 1;
end

addpath('X:\Common\Lab_Papers\spatial_sorting_intestine\data_for_paper\supporting_functions');
load('X:\Common\Lab_Papers\spatial_sorting_intestine\data_for_paper\results_7_10_2020.mat');
load('X:\Common\Lab_Papers\spatial_sorting_intestine\data_for_paper\1_SST_Protein_mRNA_TE_parsed.mat');


sst2 = sst;
sst2.protein_norm_med=sst.protein_norm_median;

yyr_mat = [];
yyp_mat = [];

[n,m] = get_subplot_size(gg);

figure;
for i = 1 : length(gg)
    clear yyr
    clear yyp
    clear yyp2
    subplot(m,n,i);
    ind =   find(strcmpi(sst.gene_name,gg{i}));
    if ~isempty(ind)
%                 yyr = sst.mRNA_norm_tans(ind,:);
        %         yyp = sst.protein_norm_median(ind,:);
        [y,t]=plot_fit_single_gene(sst2,gg{i},delta_all,score_all,0);
        yyr = smoothdata(sst.mRNA_norm_tans(ind,:),'loess',6);
        yyp = smoothdata(sst.protein_norm_median(ind,:),'loess',6);
        %         yyp2= sst.protein_norm_median(ind,:);
        yyr_mat = [yyr_mat; yyr];
        yyp_mat = [yyp_mat; yyp];
        if NORM == 1 % MEAN norm
            plot_patch(linspace(0,1,length(yyr)),yyr/mean(yyr),sst.mRNA_sem(ind,:)/mean(yyr),'b');hold on;
            plot_patch(linspace(0,1,length(yyp)),yyp/mean(yyp),sst.protein_sem(ind,:)/mean(yyp),'r'); hold on;
            plot(linspace(0,1,length(t)),y/mean(y),'k--','linewidth',2);
            %         plot_patch(linspace(0,1,length(yyp2)),yyp2/mean(yyp2),sst.protein_sem(ind,:)/mean(yyp2),'k');
            %plot(linspace(0,1,length(yyr)),cumsum(yyr)/mean(cumsum(yyr)),'g');
            %             ylabel('Relative expression [Mean norm]','FontSize',8)
        elseif NORM == 2 %MAX roem
            plot_patch(linspace(0,1,length(yyr)),yyr/max(yyr),sst.mRNA_sem(ind,:)/max(yyr),'b');
            plot_patch(linspace(0,1,length(yyp)),yyp/max(yyp),sst.protein_sem(ind,:)/max(yyp),'r');
            %         plot_patch(linspace(0,1,length(yyp2)),yyp2/max(yyp2),sst.protein_sem(ind,:)/max(yyp2),'k');
            %plot(linspace(0,1,length(yyr)),cumsum(yyr)/mean(cumsum(yyr)),'g');
            %             ylabel('Relative expression [Mean norm]','FontSize',8)
        else
            plot_patch(linspace(0,1,length(yyr)),yyr,sst.mRNA_sem(ind,:),'b');
            plot_patch(linspace(0,1,length(yyp)),yyp,sst.protein_sem(ind,:),'r');
            %         plot_patch(linspace(0,1,length(yyp2)),yyp2,sst.protein_sem(ind,:),'k');
            %plot(linspace(0,1,length(yyr)),cumsum(yyr),'g');
            %             ylabel('Relative expression','FontSize',8)
        end
        title([gg{i},', R=' num2str(score_all(ind),'%.1e') ', t_{1/2}=' num2str(log(2)./delta_all(ind),'%.1f'), ' hr']);
        xticks(linspace(0,1,length(yyr)));
        xticklabels({'V1','V2','V3','V4','V5','V6'});
%         ylim([0 max(ylim)]);
        set(gca,'FontSize',10);
        grid on;
        box on;
        num2str(sst.protein_mice_count(ind));
    end
    if i == 1
        legend('mRNA','','Protein','','dpdt','Orientation','horizontal');
    end
end