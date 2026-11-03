function Figures=study_plot_single(R,Q,F,O,output_dir)
%STUDY_PLOT_SINGLE Six figures for one selected, ten-source realization.
% All quantities come from the supplied Result. This function never refits.
% Coordinates are metres and fields are nT. Glyph sizes are schematic.
if nargin<4,O=struct();end
O=defaults(O);assert(numel(Q)==3&&numel(F.method_names)==3);
Figures=struct();
Figures.model=model_figure(R,O);
study_export_figure(Figures.model,output_dir,'Single01_source_model',O);
Figures.measured=measured_figure(R,F,O);
study_export_figure(Figures.measured,output_dir,'Single02_measured_field',O);
Figures.predicted=prediction_figure(R,Q,F,O);
study_export_figure(Figures.predicted,output_dir,'Single03_predicted_fields',O);
Figures.residuals=residual_figure(R,Q,F,O);
study_export_figure(Figures.residuals,output_dir,'Single04_residual_fields',O);
Figures.recall=recall_figure(R,Q,O);
study_export_figure(Figures.recall,output_dir,'Single05_recall_comparison',O);
Figures.positions=position_figure(R,Q,O);
study_export_figure(Figures.positions,output_dir,'Single06_source_recovery',O);
end

function O=defaults(O)
n={'visible','save_figures','save_editable_fig','font_name','dpi','model_glyph_radius_m','model_depth_exaggeration','colors'};
v={'on',true,false,'Times New Roman',300,2,5,[.36 .44 .57;.84 .55 .19;.03 .49 .55]};
for k=1:numel(n),if ~isfield(O,n{k}),O.(n{k})=v{k};end,end
if ~isfield(O,'residual_color_limit_nT'),O.residual_color_limit_nT=50;end
assert(isempty(O.residual_color_limit_nT)|| ...
    (isscalar(O.residual_color_limit_nT)&&isfinite(O.residual_color_limit_nT)&&O.residual_color_limit_nT>0));
assert(isfinite(O.model_glyph_radius_m)&&O.model_glyph_radius_m>0);
assert(isfinite(O.model_depth_exaggeration)&&O.model_depth_exaggeration>0);
end

function f=make_figure(R,O,topic,position)
f=figure('Color','w','Name',[R.Meta.name ': ' topic],'NumberTitle','off',...
    'Visible',O.visible,'Position',position,'Renderer','painters');
set(f,'DefaultAxesFontName',O.font_name,'DefaultTextFontName',O.font_name);
annotation(f,'textbox',[.04 .93 .92 .055],'String',sprintf(...
    'Selected layout %d | 10 sources | noise %g%% | %s',R.Meta.realization,100*R.Meta.noise_fraction,topic),...
    'EdgeColor','none','HorizontalAlignment','center','FontName',O.font_name,...
    'FontSize',14,'FontWeight','bold','Interpreter','none');
end

function style_axis(ax,O)
set(ax,'FontName',O.font_name,'FontSize',12,'FontWeight','bold', ...
    'TickDir','in','LineWidth',1.5,'Box','on','Layer','top');
grid(ax,'off');
end

function footnote(f,s,O)
if ~iscell(s),s={s};end
for k=1:numel(s)
    h=.07/numel(s);y=.025+(numel(s)-k)*h;
    annotation(f,'textbox',[.05 y .90 h],'String',s{k},'EdgeColor','none',...
        'HorizontalAlignment','center','VerticalAlignment','middle',...
        'FontName',O.font_name,'FontSize',10,'Interpreter','none');
end
end

function panel_letter(ax,c)
text(ax,.025,.975,['(' c ')'],'Units','normalized','FontSize',14,...
    'FontWeight','bold','VerticalAlignment','top','Color','k');
end

function f=model_figure(R,O)
f=make_figure(R,O,'true source model',[60 70 1300 650]);
P=R.Scene.Truth.position;M=R.Scene.Truth.moment;mag=sqrt(sum(M.^2,2));
roi=R.SourceInfo.source_bounds_xy;
ax=axes('Parent',f,'Position',[.055 .16 .53 .67]);hold(ax,'on');style_axis(ax,O);
% The numerical z parameter is positive depth. Physical magnetization points
% upwards when its third component is positive, hence -u_z on the depth axis.
scan=[min(R.Scene.X(:)) min(R.Scene.Y(:));max(R.Scene.X(:)) max(R.Scene.Y(:))];
sx=[scan(1,1) scan(2,1) scan(2,1) scan(1,1) scan(1,1)];
sy=[scan(1,2) scan(1,2) scan(2,2) scan(2,2) scan(1,2)];
plot3(ax,sx,sy,zeros(size(sx)),'-','Color',[.50 .50 .50],'LineWidth',1);
[xx,yy,zz]=sphere(18);r=O.model_glyph_radius_m;
for k=1:size(P,1)
    % Compensate the display-only vertical stretch so the point glyph stays round.
    surf(ax,P(k,1)+r*xx,P(k,2)+r*yy,P(k,3)+r/O.model_depth_exaggeration*zz,...
        mag(k)*ones(size(xx)),'EdgeColor','none','FaceColor','interp');
    u=M(k,:)/mag(k);a=5;
    quiver3(ax,P(k,1),P(k,2),P(k,3),a*u(1),a*u(2),-a*u(3),0,...
        'Color',[.18 .18 .18],'LineWidth',1,'MaxHeadSize',.5);
    text(ax,P(k,1)+r,P(k,2)+r,P(k,3)-.5,sprintf('%d',k),'FontSize',12,'FontWeight','bold','Color','k');
end
set(ax,'ZDir','reverse');view(ax,[-38 24]);
xlim(ax,scan(:,1).'+[-3 3]);ylim(ax,scan(:,2).'+[-3 3]);
zlim(ax,[-1 max(P(:,3))+r+2]);
daspect(ax,[1 1 1/O.model_depth_exaggeration]);
xlabel(ax,'x (m)');ylabel(ax,'y (m)');zlabel(ax,'Depth (m)');
title(ax,'True positions and magnetic-moment directions','FontWeight','normal');
colormap(ax,viridis_map());cl=finite_range(mag);caxis(ax,cl);
if exist('OCTAVE_VERSION','builtin')&&strcmp(graphics_toolkit(f),'gnuplot')
    % GNUPlot reduces the native 3-D plotting rectangle a second time after
    % applying a manual data aspect ratio. Enlarge only its drawing viewport;
    % keep physical limits and the requested depth aspect ratio unchanged.
    set(ax,'Position',[-.105 -.005 .85 1.0]);title(ax,'');
    annotation(f,'textbox',[.01 .83 .58 .055],...
        'String','True positions and magnetic-moment directions',...
        'EdgeColor','none','HorizontalAlignment','center','FontSize',12,...
        'FontName',O.font_name);
end
% Second panel gives an undistorted, easy-to-read planar reference.
bx=axes('Parent',f,'Position',[.64 .22 .26 .53]);hold(bx,'on');style_axis(bx,O);
scatter(bx,P(:,1),P(:,2),180,mag,'filled');
for k=1:size(P,1),text(bx,P(k,1)+1.5,P(k,2)+2.0,sprintf('%d',k),'FontSize',12,'FontWeight','bold');end
rectangle('Parent',bx,'Position',[roi(1,:) diff(roi,1,1)],'EdgeColor',[.55 .55 .55],'LineStyle','--');
axis(bx,'equal');xlim(bx,roi(:,1).'+[-5 5]);ylim(bx,roi(:,2).'+[-5 5]);
xlabel(bx,'x (m)');ylabel(bx,'y (m)');title(bx,'Plan view','FontWeight','normal');
colormap(bx,viridis_map());caxis(bx,cl);p=get(bx,'Position');cb=colorbar(bx);set(bx,'Position',p);set(cb,'Position',[.915 .22 .012 .53]);
ylabel(cb,'Moment magnitude (A m^2)');
footnote(f,sprintf('Point-dipole forward model; sphere sizes are schematic. Depth display exaggeration: %gx. The outlined plane is the measurement surface.',...
    O.model_depth_exaggeration),O);
end

function f=measured_figure(R,F,O)
f=make_figure(R,O,'clean and measured magnetic fields',[70 70 1220 630]);
cl=finite_range([F.clean_nT(:);F.measured_nT(:)]);
for j=1:2
    ax=axes('Parent',f,'Position',[.075+(j-1)*.46 .19 .34 .65]);
    if j==1,z=F.clean_nT;ttl='Noise-free B_z';else,z=F.measured_nT;ttl='Measured B_z';end
    draw_field(ax,F,z,cl,jet(256),O);title(ax,ttl,'FontWeight','normal');panel_letter(ax,char('a'+j-1));
    if j==2,p=get(ax,'Position');cb=colorbar(ax);set(ax,'Position',p);set(cb,'Position',[.91 .19 .015 .65]);ylabel(cb,'B_z (nT)');end
end
footnote(f,sprintf('Same color limits in both panels. Gaussian noise scale = %.4g nT (%g%% of the clean-map RMS).',...
    R.Scene.Sigma_nT,100*R.Meta.noise_fraction),O);
end

function f=prediction_figure(R,Q,F,O)
f=make_figure(R,O,'three-method predicted magnetic fields',[80 80 1400 550]);
v=F.measured_nT(:);for j=1:3,v=[v;F.prediction_nT{j}(:)];end %#ok<AGROW>
cl=finite_range(v);
for j=1:3
    ax=axes('Parent',f,'Position',[.055+(j-1)*.30 .21 .245 .62]);
    draw_field(ax,F,F.prediction_nT{j},cl,jet(256),O);
    title(ax,Q(j).method_label,'FontWeight','normal');panel_letter(ax,char('a'+j-1));
    if j==3,p=get(ax,'Position');cb=colorbar(ax);set(ax,'Position',p);set(cb,'Position',[.935 .21 .012 .62]);ylabel(cb,'Predicted B_z (nT)');end
end
footnote(f,'Common color limits. Window-method predictions include their saved fitted affine background; no additional fitting is performed for plotting.',O);
end

function f=residual_figure(R,Q,F,O)
f=make_figure(R,O,'observed minus predicted magnetic fields',[90 80 1400 580]);
a=0;for j=1:3,a=max(a,max(abs(F.residual_nT{j}(:))));end;a=max(a,eps);
if ~isempty(O.residual_color_limit_nT),a=O.residual_color_limit_nT;end
for j=1:3
    ax=axes('Parent',f,'Position',[.055+(j-1)*.30 .23 .245 .59]);
    draw_field(ax,F,F.residual_nT{j},[-a a],residual_map(),O);
    title(ax,sprintf('%s | RMS %.4g nT',Q(j).method_label,Q(j).residual_RMSE_nT),'FontWeight','normal','FontSize',10);
    panel_letter(ax,char('a'+j-1));
    if j==3,p=get(ax,'Position');cb=colorbar(ax);set(ax,'Position',p);set(cb,'Position',[.935 .23 .012 .59]);ylabel(cb,'Measured - predicted (nT)');end
end
footnote(f,sprintf('Common display limits: +/- %g nT; values outside the range are color-saturated. RMS uses all unclipped residual values.',a),O);
end

function draw_field(ax,F,z,cl,cmap,O)
imagesc(ax,F.X(1,:),F.Y(:,1),z);set(ax,'YDir','normal');axis(ax,'image');
style_axis(ax,O);grid(ax,'off');colormap(ax,cmap);caxis(ax,cl);
xlabel(ax,'x (m)');ylabel(ax,'y (m)');
end

function f=recall_figure(R,Q,O)
f=make_figure(R,O,'recall and source-set evaluation',[90 60 1250 630]);
ax=axes('Parent',f,'Position',[.075 .25 .49 .58]);hold(ax,'on');style_axis(ax,O);
for j=1:3
    val=100*Q(j).Recall;
    if val>0,rectangle('Parent',ax,'Position',[j-.27 0 .54 val],...
            'FaceColor',O.colors(j,:),'EdgeColor',O.colors(j,:));end
    text(ax,j,val+3,sprintf('%g%%  (%d/10)',val,Q(j).TP),...
        'HorizontalAlignment','center','FontSize',11,'FontWeight','bold');
end
set(ax,'XTick',1:3,'XTickLabel',{'Primin-based','Ind. + fusion','Window game'});
xlim(ax,[.4 3.6]);ylim(ax,[0 114]);ylabel(ax,'Recall (%)');
title(ax,'Matched true sources / 10','FontWeight','normal');
bx=axes('Parent',f,'Position',[.63 .24 .33 .58],'Visible','off');hold(bx,'on');
xlim(bx,[0 1]);ylim(bx,[0 1]);
heads={'','Primin','Fusion','Game'};xs=[.01 .39 .64 .88];
for j=1:4,text(bx,xs(j),.94,heads{j},'FontWeight','bold','HorizontalAlignment','center','FontSize',11);end
labs={'TP','FP','FN','F1','GOSPA (m)'};
V=[[Q.TP];[Q.FP];[Q.FN];[Q.F1];[Q.GOSPA_m]];
for k=1:5
    y=.81-(k-1)*.15;text(bx,.01,y,labs{k},'HorizontalAlignment','left','FontSize',11);
    for j=1:3
        if k<=3,s=sprintf('%d',V(k,j));else,s=sprintf('%.3f',V(k,j));end
        text(bx,xs(j+1),y,s,'HorizontalAlignment','center','FontSize',11,'Color',O.colors(j,:));
    end
end
plot(bx,[0 1],[.88 .88],'-','Color',[.5 .5 .5]);
footnote(f,sprintf('One-to-one 3D matching within %g m. Standard GOSPA: p = 1, c = %g m, alpha = 2. Single-case values, without uncertainty bars.',...
    Q(1).position_tolerance_m,Q(1).position_tolerance_m),O);
end

function f=position_figure(R,Q,O)
f=make_figure(R,O,'source recovery',[60 80 1400 580]);
P=R.Scene.Truth.position;roi=R.SourceInfo.source_bounds_xy;allxy=[P(:,1:2);roi];
for j=1:3,M=get_method(R,Q(j).method_name);allxy=[allxy;M.Fit.position(:,1:2)];end %#ok<AGROW>
lo=min(allxy,[],1)-5;hi=max(allxy,[],1)+5;width=max(hi-lo);center=(lo+hi)/2;limits=[center-width/2;center+width/2];
c=Q(1).position_tolerance_m;lastax=[];markers={'o','d','s'};
for j=1:3
    ax=axes('Parent',f,'Position',[.055+(j-1)*.295 .23 .245 .60]);hold(ax,'on');axis(ax,'equal');style_axis(ax,O);
    xlim(ax,limits(:,1).');ylim(ax,limits(:,2).');xlabel(ax,'x (m)');ylabel(ax,'y (m)');title(ax,Q(j).method_label,'FontWeight','normal');
    rectangle('Parent',ax,'Position',[roi(1,:) diff(roi,1,1)],'EdgeColor',[.60 .60 .60],'LineStyle','--');
    M=get_method(R,Q(j).method_name);E=M.Fit.position;D=M.Details;ti=D.matched_truth_indices(:);ei=D.matched_estimate_indices(:);
    for k=1:numel(ti)
        plot(ax,[P(ti(k),1) E(ei(k),1)],[P(ti(k),2) E(ei(k),2)],'-','Color',[.70 .70 .70],'LineWidth',.7);
        cmap=viridis_map();ci=1+round(255*min(1,max(0,D.matched_errors(k)/c)));
        plot(ax,E(ei(k),1),E(ei(k),2),markers{j},'MarkerSize',8,'MarkerFaceColor',cmap(ci,:),...
            'MarkerEdgeColor',[.2 .2 .2],'LineWidth',.7);
    end
    ui=D.unmatched_estimate_indices(:);if ~isempty(ui),plot(ax,E(ui,1),E(ui,2),'+','MarkerSize',10,'LineWidth',1.7,'Color',[.45 .45 .45]);end
    plot(ax,P(:,1),P(:,2),'x','Color',[.85 .12 .12],'MarkerSize',8,'LineWidth',1.8);
    for k=1:size(P,1),text(ax,P(k,1)+1,P(k,2)+2,sprintf('%d',k),'FontSize',12,'FontWeight','bold','Color','k');end
    panel_letter(ax,char('a'+j-1));
    text(ax,.5,.04,sprintf('TP %d   FP %d   FN %d',Q(j).TP,Q(j).FP,Q(j).FN),'Units','normalized',...
        'HorizontalAlignment','center','FontSize',10,'BackgroundColor','w');
    colormap(ax,viridis_map());caxis(ax,[0 c]);lastax=ax;
end
p=get(lastax,'Position');cb=colorbar(lastax);set(lastax,'Position',p);set(cb,'Position',[.93 .25 .013 .54]);
ylabel(cb,'3D position error (m)');
footnote(f,{'Red x: true source; colored markers: matched estimates; gray +: unmatched estimates; gray lines: matched pairs.',...
    'All estimates are shown. Colors encode three-dimensional matching error, although positions are displayed in plan view.'},O);
end

function M=get_method(R,name)
M=[];A=R.Methods;if ~iscell(A),A=num2cell(A);end
for k=1:numel(A),if strcmp(A{k}.name,name),M=A{k};return;end,end
end

function cl=finite_range(v)
v=v(isfinite(v));if isempty(v),cl=[0 1];else,cl=[min(v) max(v)];end
if cl(2)<=cl(1),cl=cl+[-1 1]*max(1,abs(cl(1))*.01);end
end

function c=viridis_map()
a=[.267 .005 .329;.254 .265 .530;.164 .471 .558;.135 .659 .518;.478 .821 .318;.993 .906 .144];
c=interp1(linspace(0,1,size(a,1)),a,linspace(0,1,256),'linear');
end

function c=residual_map()
c=interp1([0 .5 1],[.14 .31 .69;1 1 1;.75 .10 .15],linspace(0,1,256),'linear');
end
