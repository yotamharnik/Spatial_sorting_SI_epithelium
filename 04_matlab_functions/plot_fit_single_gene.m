function [y,t]=plot_fit_single_gene(sst,gene_name,delta_vec,score_vec,DISP)

if nargin<5
    DISP=1;
end
ind=find(strcmpi(sst.gene_name,gene_name));
if isempty(ind)
    return;
end
delta=delta_vec(ind);
tspan=0:72;

m=sst.mRNA_norm(ind,:);
p=sst.protein_norm(ind,:);
m=m/mean(m);
p=p/mean(p);
p2=sst.protein_norm_med(ind,:);
p2=p2/mean(p2);

xx=(0:72)';
m=interp1(linspace(0,72,6),m,xx);
p=interp1(linspace(0,72,6),p,xx);
p2=interp1(linspace(0,72,6),p2,xx);

mt=[xx m];
save mdat mt delta

[t,y] = ode23(@mRNA2prot,tspan,p(1));

if DISP
    subplot(1,2,1);
    plot(mt(:,1),mt(:,2)/mean(mt(:,2)),'linewidth',2);
    hold on;
    plot(t,y/mean(y),'r','linewidth',2)
    plot(t,p/mean(p),'m','linewidth',2)
    plot(t,p2/mean(p2),'k','linewidth',2)
    
    
    legend('mRNA',['\delta=' num2str(delta)],'protein');
    box on;
    xlabel('Time (hrs)');
    ylabel('Mean-normalized expression');
    title([gene_name ', R=' num2str(score_vec(ind)) ', t_{1/2}=' num2str(log(2)./delta,'%.1f') ' hr']);
    
    subplot(1,2,2);
    m=sst.mRNA_norm(ind,:);
    p=sst.protein_norm(ind,:);
    p2=sst.protein_norm_med(ind,:);
    
    mse=sst.mRNA_sem(ind,:);
    pse=sst.protein_sem(ind,:);
    plot_patch(linspace(0,1,length(m)),m,mse,'b');
    hold on;
    xlabel('zone');
    ylabel('Expression');
    
    plot(linspace(0,1,length(y)),y*mean(p)./mean(y),'r','linewidth',2)
    plot_patch(linspace(0,1,length(p)),p,pse,'m')
    plot(linspace(0,1,length(p2)),p2,'m')
    
    box on;
    
    set(gcf,'position',[159         318        1261         480]);
end