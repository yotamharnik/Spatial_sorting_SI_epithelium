function [yyr_mat,yyp_mat] = plot_sst_zonation_m1_m4(gg,NORM,plot_flag);

if nargin <2
    NORM = 0;
    plot_flag = 1;
end
if nargin <3
    plot_flag = 1;
end

load('X:\Yotam\matlab_projects\spatial_sorting_thesis\0_data\1_tidy_data\6_SST_Protein_mRNA_TE_parsed_ver2.mat');
load('X:\Yotam\matlab_projects\spatial_sorting_thesis\0_data\1_tidy_data\3_transcription_efficiency_vs_protein_stability_m1_m4.mat');
addpath('X:\Yotam\matlab_projects\spatial_sorting_thesis\1_code\2_functions');

yyr_mat = [];
yyp_mat = [];

gg = intersect(lower(sst.gene_name),lower(gg));

[n,m] = get_subplot_size(gg);

if plot_flag
    figure;
    suptitle('blue - mRNA, red - Protein');
end
for i = 1 : length(gg)
    clear yyr
    clear yyp
    clear yyp2
    subplot(n,m,i);
    ind =   find(strcmpi(sst.gene_name,gg{i}));
    indTE = find(strcmpi(TE_PS.gene_name,gg{i}));
    if ~isempty(ind)
%         yyr = sst.mRNA_norm(ind,:);
%         yyp = sst.protein_norm_median(ind,:);
        yyr = smoothdata(sst.mRNA_norm_tans(ind,:),'loess',5);
        yyp = smoothdata(sst.protein_norm_median(ind,:),'loess',5);
%         yyp2= sst.protein_norm_median(ind,:);
        yyr_mat = [yyr_mat; yyr];
        yyp_mat = [yyp_mat; yyp];
    if NORM == 1 % MEAN norm
        plot_patch(linspace(0,1,length(yyr)),yyr/mean(yyr),sst.mRNA_sem(ind,:)/mean(yyr),'b');
        plot_patch(linspace(0,1,length(yyp)),yyp/mean(yyp),sst.protein_sem(ind,:)/mean(yyp),'r');
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
    if ~isempty(indTE)
        hold on;
        %text(0.2,0.2,);
    end
    if ~isempty(indTE)
%         title([gg{i}, ' n(p)= ',num2str(sst.protein_mice_count(ind)) ', TE%=' num2str(100*length(find(TE_PS.mat(:,1)>TE_PS.mat(indTE,1)))/size(TE_PS.mat,1))]);
        title(gg{i},'FontSize',12);
    else
        title(gg{i},'FontSize',12);
%         title([gg{i},' n(p)= ',num2str(sst.protein_mice_count(ind))]);
    end
    %         xlabel('Zone');
    xticks(linspace(0,1,length(yyr)));
    xticklabels({'V1','V2','V3','V4','V5','V6'});
    ylim([0 max(ylim)]);
    set(gca,'FontSize',10);
    grid on;
    box on;
    num2str(sst.protein_mice_count(ind));
end
end