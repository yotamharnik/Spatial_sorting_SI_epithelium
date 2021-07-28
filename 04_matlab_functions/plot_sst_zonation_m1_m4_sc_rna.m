function [yyr_mat,yyp_mat] = plot_sst_zonation_m1_m4_sc_rna(gg,NORM,plot_flag);

if nargin == 1
    NORM = 0;
    plot_flag = 1;
end
if nargin == 2
    plot_flag = 1;
end

load('X:\Common\Lab_Papers\spatial_sorting_intestine\data_for_paper\2_sst_with_single_cell.mat');
addpath('X:\Yotam\matlab_projects\spatial_sorting_thesis\1_code\2_functions');

yyr_mat = [];
yyp_mat = [];

gg = intersect(lower(gg),lower(sst.gene_name));

[n,m] = get_subplot_size(gg);


for i = 1 : length(gg)
    clear yyr
    clear yyp
    clear yyp2
    subplot(n,m,i);
    ind =   find(strcmpi(sst.gene_name,gg{i}));
    if ~isempty(ind)
        %         yyr = sst.sc_mn(ind,:);
        %         yyp = sst.protein_norm_median(ind,:);
%         yyr = smoothdata(sst.sc_mn(ind,:),'loess',6);
%         yyp = smoothdata(sst.protein_norm(ind,:),'loess',6);
        yyr = sst.sc_mn(ind,:);
        yyp = sst.protein_norm(ind,:);
        %         yyp2= sst.protein_norm_median(ind,:);
        yyr_mat = [yyr_mat; yyr];
        yyp_mat = [yyp_mat; yyp];
        if plot_flag
            if NORM == 1 % MEAN norm
                plot_patch(linspace(0,1,length(yyr)),yyr/mean(yyr),sst.sc_se(ind,:)/mean(yyr),'b');
                plot_patch(linspace(0,1,length(yyp)),yyp/mean(yyp),sst.protein_sem(ind,:)/mean(yyp),'r');
                %         plot_patch(linspace(0,1,length(yyp2)),yyp2/mean(yyp2),sst.protein_sem(ind,:)/mean(yyp2),'k');
                %plot(linspace(0,1,length(yyr)),cumsum(yyr)/mean(cumsum(yyr)),'g');
                %             ylabel('Relative expression [Mean norm]','FontSize',8)
            elseif NORM == 2 %MAX roem
                plot_patch(linspace(0,1,length(yyr)),yyr/max(yyr),sst.sc_se(ind,:)/max(yyr),'b');
                plot_patch(linspace(0,1,length(yyp)),yyp/max(yyp),sst.protein_sem(ind,:)/max(yyp),'r');
                %         plot_patch(linspace(0,1,length(yyp2)),yyp2/max(yyp2),sst.protein_sem(ind,:)/max(yyp2),'k');
                %plot(linspace(0,1,length(yyr)),cumsum(yyr)/mean(cumsum(yyr)),'g');
                %             ylabel('Relative expression [Mean norm]','FontSize',8)
            else
                plot_patch(linspace(0,1,length(yyr)),yyr,sst.sc_se(ind,:),'b');
                plot_patch(linspace(0,1,length(yyp)),yyp,sst.protein_sem(ind,:),'r');
                %         plot_patch(linspace(0,1,length(yyp2)),yyp2,sst.protein_sem(ind,:),'k');
                %plot(linspace(0,1,length(yyr)),cumsum(yyr),'g');
                %             ylabel('Relative expression','FontSize',8)
            end
            title(gg{i},'FontSize',12);
            xticks(linspace(0,1,length(yyr)));
            xticklabels({'V1','V2','V3','V4','V5','V6'});
            ylim([0 max(ylim)]);
            set(gca,'FontSize',10);
            grid on;
            box on;
            num2str(sst.protein_mice_count(ind));
        else
            continue
        end
    end
end
end