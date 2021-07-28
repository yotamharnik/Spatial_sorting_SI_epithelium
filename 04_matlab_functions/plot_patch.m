function plot_patch(x,mn,se,color,PATCH_ONLY)

% Plots a patch with mean+- one standard error and the given color
% PTCH_ONLY=0 - draw both, 1 - draw only patch, 2 - draw only line 
% x=0:20;
% y1=smooth(mn1(21:end));
% y2=smooth(se1(21:end));

mn(isnan(mn))=0;
se(isnan(se))=0;
if nargin<5,
    PATCH_ONLY=0;
end
x=x(:);mn=mn(:);se=se(:);
hold on;
if PATCH_ONLY==0,
    plot(x,mn,'LineWidth',2,'color',color);
    patch([x ;x(end:-1:1)],[mn+se; mn(end:-1:1)-se(end:-1:1)],color,'EdgeColor','none','FaceAlpha',0.1);
elseif PATCH_ONLY==1,
    patch([x ;x(end:-1:1)],[mn+se; mn(end:-1:1)-se(end:-1:1)],color,'EdgeColor','none','FaceAlpha',0.1);
else
    plot(x,mn,'LineWidth',3,'color',color);
end
%axis square;
%axis([xlim 0 1.1]);
set(gca,'FontSize',16);