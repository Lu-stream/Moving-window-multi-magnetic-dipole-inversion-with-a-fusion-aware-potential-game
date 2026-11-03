function Figures = plot_area1_paper(Results, output_dir, Plot)
%PLOT_AREA1_PAPER Paper figures from saved measured results; no inversion.
% MATLAB R2021b compatible. CT symbols are geometric references, not known
% magnetic dipole centers. All matching uses the saved 3-D evaluation.
% Field maps are displayed in microtesla (saved maps are in nanotesla).
% Matched estimates use the method color; unmatched estimates are gray +.
% No colorbar is drawn for the position plot because color is categorical.
if nargin < 3, Plot=struct(); end
Plot=defaults(Plot);
if ~exist(output_dir,'dir'), mkdir(output_dir); end
assert(isfield(Results,'Methods') && numel(Results.Methods)==3,...
 'Expected the three measured-data methods in Results.Methods.');
offset=Results.Metadata.xy_offset_um(:)';
X=Results.Data.X+offset(1); Y=Results.Data.Y+offset(2);
B=Results.Data.Bz/1000;
extent=[min(X(:)) max(X(:)) min(Y(:)) max(Y(:))];
refs=Results.Reference.position; refs(:,1:2)=bsxfun(@plus,refs(:,1:2),offset);
ids=Results.Reference.ids(:);
labels={'Primin-based','Independent + fusion','Fusion-aware game'};
colors=[.36 .44 .57; .84 .55 .19; .03 .49 .55];
markers={'o','d','s'};
prediction=cell(1,3); residual=cell(1,3); residual_max=0;
for j=1:3
 M=Results.Methods{j};
 prediction{j}=M.PhysicalInfo.prediction/1000;
 residual{j}=M.PhysicalInfo.residual/1000;
 assert(isequal(size(residual{j}),size(B)), 'Residual grid size mismatch.');
 err=B-prediction{j}-residual{j};
 assert(max(abs(err(:)))<1e-8, 'Saved residual is not measured minus predicted.');
 residual_max=max(residual_max,max(abs(residual{j}(:))));
end
if isempty(Plot.residual_limit_uT), residual_limit=max(residual_max,eps);
else
 residual_limit=Plot.residual_limit_uT;
 assert(isscalar(residual_limit)&&isfinite(residual_limit)&&residual_limit>0,...
  'Plot.residual_limit_uT must be empty or a positive scalar.');
end
% Keep <=256 entries for Octave/gnuplot's indexed-image PDF export.
cmap=diverging_map(256);
Figures=struct();

%% Registered CT geometry and measured Bz.
Figures.input=make_figure('area1_paper_input',[80 80 1130 410],Plot);
ax=axes('Parent',Figures.input,'Position',[.065 .16 .375 .76]); hold(ax,'on');
cuboid_file=fullfile(Plot.input_dir,'FWInput-FineCuboids-A1.txt');
if exist(cuboid_file,'file')==2
 C=dlmread(cuboid_file);
 % Columns 4 and 5 are half widths, as documented by the source archive.
 % Draw each grain's cuboids as one patch collection to keep rendering light.
 grain_colors=lines(numel(ids));
 for k=1:numel(ids)
  q=C(C(:,7)==ids(k),:);
  xx=[q(:,1)-q(:,4),q(:,1)+q(:,4),q(:,1)+q(:,4),q(:,1)-q(:,4)]';
  yy=[q(:,2)-q(:,5),q(:,2)-q(:,5),q(:,2)+q(:,5),q(:,2)+q(:,5)]';
  patch('Parent',ax,'XData',xx,'YData',yy,'FaceColor',grain_colors(k,:),...
   'EdgeColor','none');
 end
end
plot_refs(ax,refs,ids,Plot);
style_axes(ax,extent,Plot);
panel_label(ax,'(a)',Plot);
if Plot.show_titles, title(ax,'CT grain geometry','FontWeight','bold'); end
input_map_pos=[.515 .16 .375 .76];
ax=axes('Parent',Figures.input,'Position',input_map_pos);
imagesc(ax,X(1,:),Y(:,1),B); set(ax,'YDir','normal');
colormap(ax,cmap); field_limit=max(max(abs(B(:))),eps); caxis(ax,[-field_limit field_limit]);
style_axes(ax,extent,Plot); hold(ax,'on'); panel_label(ax,'(b)',Plot);
cb=colorbar(ax); set(ax,'Position',input_map_pos);
set(cb,'Position',[.910 .245 .014 .57]);
style_colorbar(cb,'',Plot);
method_label(ax,'Measured B_z (\muT)',Plot,1.02);
if Plot.show_titles, title(ax,'Measured field','FontWeight','bold'); end
save_one(Figures.input,output_dir,'Fig_measured_input',Plot);

%% Full-map residuals: identical limits, metrics from unclipped residuals.
Figures.residual=make_figure('area1_paper_residual',[100 100 1440 440],Plot);
for j=1:3
 residual_pos=[.055+(j-1)*.300 .23 .24 .63];
 ax=axes('Parent',Figures.residual,'Position',residual_pos);
 imagesc(ax,X(1,:),Y(:,1),residual{j}); set(ax,'YDir','normal');
 colormap(ax,cmap); caxis(ax,[-residual_limit residual_limit]);
 style_axes(ax,extent,Plot); hold(ax,'on');
 if j>1,ylabel(ax,'');end
 panel_label(ax,sprintf('(%c)',char(96+j)),Plot);
 method_label(ax,labels{j},Plot);
 rms_uT=Results.Methods{j}.PhysicalInfo.RMSE/1000;
 text(ax,.5,-.24,sprintf('RMS = %.3f \\muT',rms_uT),...
  'Units','normalized','HorizontalAlignment','center',...
  'FontName',Plot.font_name,'FontSize',Plot.axis_font_size,...
  'FontWeight','bold','Interpreter','tex');
end
cb=colorbar(ax); set(ax,'Position',residual_pos);
set(cb,'Position',[.908 .28 .014 .49]);
style_colorbar(cb,'Measured - predicted B_z (\muT)',Plot);
save_one(Figures.residual,output_dir,'Fig_measured_residuals',Plot);

%% Localization: retain every output and the fixed 3-D assignment.
Figures.positions=make_figure('area1_paper_positions',[110 110 1440 480],Plot);
for j=1:3
 ax=axes('Parent',Figures.positions,'Position',[.055+(j-1)*.322 .29 .275 .59]);
 hold(ax,'on'); M=Results.Methods{j};
 est=M.Fit.position; est(:,1:2)=bsxfun(@plus,est(:,1:2),offset);
 matched=M.Details.matched_estimate_indices(:);
 unmatched=M.Details.unmatched_estimate_indices(:);
 if ~isempty(unmatched)
  plot(ax,est(unmatched,1),est(unmatched,2),'+','Color',[.42 .42 .42],...
   'LineStyle','none','MarkerSize',9,'LineWidth',1.8);
 end
 if ~isempty(matched)
  plot(ax,est(matched,1),est(matched,2),markers{j},'Color',[.12 .12 .12],...
   'MarkerFaceColor',colors(j,:),'LineStyle','none','MarkerSize',9,'LineWidth',1.2);
 end
 plot_refs(ax,refs,ids,Plot); style_axes(ax,extent,Plot);
 if j>1,ylabel(ax,'');end
 panel_label(ax,sprintf('(%c)',char(96+j)),Plot); method_label(ax,labels{j},Plot);
 text(ax,.5,-.28,sprintf('Matched %d   Unmatched %d   Missed %d',...
  M.Metrics.TP,M.Metrics.FP,M.Metrics.FN),...
  'Units','normalized','HorizontalAlignment','center','FontName',Plot.font_name,...
  'FontSize',max(10,Plot.axis_font_size-1),'FontWeight','bold');
end
annotation(Figures.positions,'textbox',[.08 .015 .84 .055],...
 'String','Red x: CT center     Filled symbol: matched estimate     Gray +: unmatched estimate',...
 'HorizontalAlignment','center','VerticalAlignment','middle','EdgeColor','none',...
 'FontName',Plot.font_name,'FontSize',max(10,Plot.axis_font_size-1),'FontWeight','bold');
save_one(Figures.positions,output_dir,'Fig_measured_localization',Plot);

if Plot.make_counts
 Figures.counts=make_figure('area1_paper_counts',[130 130 700 480],Plot);
 ax=axes('Parent',Figures.counts,'Position',[.13 .20 .82 .72]);
 counts=zeros(3,2);
 for j=1:3
  counts(j,:)=[Results.Methods{j}.Metrics.TP Results.Methods{j}.Metrics.FP];
 end
 bars=bar(ax,1:3,counts,'stacked');
 set(bars(1),'FaceColor',[.03 .49 .55]); set(bars(2),'FaceColor',[.72 .72 .72]);
 set(ax,'XTick',1:3,'XTickLabel',labels,'FontName',Plot.font_name,...
  'FontSize',Plot.axis_font_size,'FontWeight','bold','LineWidth',Plot.border_width,...
  'Box','on','TickDir','in');
 ylabel(ax,'Number of outputs','FontWeight','bold','FontSize',Plot.label_font_size);
 legend(ax,{'Matched CT reference','Unmatched output'},'Location','northwest');
 save_one(Figures.counts,output_dir,'Fig_measured_output_counts',Plot);
end

fprintf('Plotting only: %d CT centers, 3 methods, matching tolerance %.6g um.\n',...
 numel(ids),Results.Metadata.evaluation_tolerance_um);
fprintf('Residual color limits: [%.6g, %.6g] uT; full-range maximum: %.6g uT.\n',...
 -residual_limit,residual_limit,residual_max);
if residual_limit<residual_max
 fprintf('Display is saturated; the printed RMS is computed from full residuals.\n');
end
end

function P=defaults(P)
D=struct('visible','on','show_titles',false,'save_figures',false,...
 'export_dpi',600,'font_name','Times New Roman','axis_font_size',14,...
 'label_font_size',16,'border_width',2.5,'residual_limit_uT',[],...
 'input_dir','','make_counts',false,'save_editable_fig',false);
f=fieldnames(D);
for k=1:numel(f), if ~isfield(P,f{k}),P.(f{k})=D.(f{k});end,end
end

function fig=make_figure(tag,pos,P)
% Close only this script's earlier figure of the same kind, not user figures.
old=findall(0,'Type','figure','Tag',tag);
if ~isempty(old),close(old);end
fig=figure('Color','w','Name',tag,'NumberTitle','off','Tag',tag,...
 'Units','pixels','Position',pos,'Visible',P.visible,'Renderer','painters');
paper_size=pos(3:4)/100;
set(fig,'PaperUnits','inches','PaperSize',paper_size,...
 'PaperPosition',[0 0 paper_size],'PaperPositionMode','manual');
end

function style_axes(ax,extent,P)
axis(ax,extent); axis(ax,'equal'); axis(ax,extent);
set(ax,'Box','on','Layer','top','FontName',P.font_name,...
 'FontSize',P.axis_font_size,'FontWeight','bold','LineWidth',P.border_width,...
 'TickDir','in','XGrid','off','YGrid','off');
xlabel(ax,'x (\mum)','FontName',P.font_name,'FontSize',P.label_font_size,'FontWeight','bold');
ylabel(ax,'y (\mum)','FontName',P.font_name,'FontSize',P.label_font_size,'FontWeight','bold');
end

function plot_refs(ax,refs,ids,P)
plot(ax,refs(:,1),refs(:,2),'x','Color',[.88 .12 .08],...
 'LineStyle','none','MarkerSize',9,'LineWidth',2.2);
% Modest fixed text offsets move labels only; no source coordinates move.
dx=ones(size(ids))*5; dy=ones(size(ids))*5;
dy(ids==6)=-12; dx(ids==6)=-9; dy(ids==7)=7; dx(ids==7)=6;
for k=1:numel(ids)
 text(ax,refs(k,1)+dx(k),refs(k,2)+dy(k),sprintf('%d',ids(k)),...
  'FontName',P.font_name,'FontSize',P.axis_font_size,'FontWeight','bold',...
  'Color',[.1 .1 .1],'BackgroundColor','w','Margin',.1,'Clipping','on');
end
end

function panel_label(ax,s,P)
text(ax,.025,.945,s,'Units','normalized','FontName',P.font_name,...
 'FontSize',P.label_font_size,'FontWeight','bold','Color','k',...
 'BackgroundColor','w','Margin',.5,'VerticalAlignment','top');
end

function method_label(ax,s,P,y)
if nargin<4,y=1.075;end
text(ax,.5,y,s,'Units','normalized','FontName',P.font_name,...
 'FontSize',P.axis_font_size,'FontWeight','bold','HorizontalAlignment','center');
end

function style_colorbar(cb,s,P)
set(cb,'FontName',P.font_name,'FontSize',P.axis_font_size,...
 'FontWeight','bold','LineWidth',P.border_width,'TickDir','in');
if exist('OCTAVE_VERSION','builtin')
 h=ylabel(cb,s,'FontName',P.font_name,'FontSize',P.label_font_size,'FontWeight','bold');
else
 h=cb.Label;
 set(h,'String',s,'FontName',P.font_name,'FontSize',P.label_font_size,'FontWeight','bold');
end
% Explicit label placement avoids the wide auto-padding of a narrow bar.
set(h,'Units','normalized','Position',[4.2 .5 0]);
end

function map=diverging_map(n)
anchors=[.13 .27 .55; .55 .73 .87; 1 1 1; .90 .59 .48; .63 .08 .14];
map=interp1(linspace(0,1,size(anchors,1)),anchors,linspace(0,1,n));
end

function save_one(fig,out,name,P)
if ~P.save_figures,return;end
% Export is opt-in: opening the main script does not render 600-dpi images.
drawnow;
if exist('OCTAVE_VERSION','builtin')
 print(fig,fullfile(out,[name '.pdf']),'-dpdf','-painters');
 print(fig,fullfile(out,[name '.png']),'-dpng',sprintf('-r%d',P.export_dpi));
else
 exportgraphics(fig,fullfile(out,[name '.pdf']),'ContentType','vector','BackgroundColor','white');
 exportgraphics(fig,fullfile(out,[name '.png']),'Resolution',P.export_dpi,'BackgroundColor','white');
end
if P.save_editable_fig && ~exist('OCTAVE_VERSION','builtin')
 savefig(fig,fullfile(out,[name '.fig']));
end
end
