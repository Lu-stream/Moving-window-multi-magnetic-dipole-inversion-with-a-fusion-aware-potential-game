function Figures=study_plot_ensemble(Data,Summary,Options,output_dir)
%STUDY_PLOT_ENSEMBLE Selected-layout boxplots and paired noise-level summaries.
% Box quartiles use linear interpolation h=1+(n-1)*p (Hyndman-Fan type 7).
% Whiskers reach observed points within Q1-1.5*IQR and Q3+1.5*IQR.
% Every layout is drawn as a point; the three median centers are connected.
% Noise curves are per-level means, with optional +/- one sample SD,
% not confidence intervals. No data are pooled across noise levels.
if nargin<3||isempty(Options),Options=struct();end
if nargin<4,output_dir='';end
Options=defaults(Options);
assert(numel(Data.method_names)==3,'Expected the three specified methods.');
nl=numel(Data.layout_ids);nn=numel(Data.noise_fractions);
ib=find(abs(Data.noise_fractions-Options.box_noise_fraction)<=1e-12,1);
assert(~isempty(ib),'box_noise_fraction must be among the computed noise levels.');
assert(nl>=2,'At least two layouts are needed for a distribution plot.');
colors=Options.colors;labels=Data.method_labels;
Figures=struct();

% Three standalone distributions at one explicitly chosen noise level.
% FP remains in the numerical tables and the supporting noise curves below.
fields={'Recall','F1','GOSPA_m'};scales=[100 100 1];
titles={'Recall','F1 score','Standard GOSPA'};
ylabels={'Recall (%)','F_1 (%)','GOSPA (m)'};
figure_fields={'box_recall','box_f1','box_gospa'};
file_names={'Ensemble01a_recall_boxplot','Ensemble01b_f1_boxplot','Ensemble01c_gospa_boxplot'};
for k=1:3
    fig=make_figure(titles{k},[100+30*k 100 900 650],Options);
    ax=axes('Parent',fig,'Position',[.12 .23 .83 .68]);hold(ax,'on');
    values=reshape(Data.(fields{k})(:,ib,:),nl,3)*scales(k);
    draw_paired_boxes(ax,values,colors,Options.pair_lines);
    set(ax,'XTick',1:3,'XTickLabel',labels,'XLim',[.5 3.5]);
    ylabel(ax,ylabels{k});
    if Options.show_titles,title(ax,titles{k},'FontWeight','bold');end
    if k<=2
        set(ax,'YLim',[-3 106],'YTick',0:20:100);
    else
        ymax=max(values(:));set(ax,'YLim',[0 max(1,ymax)*1.13]);
    end
    style_axes(ax,Options);
    note(fig,{sprintf('%d selected layouts; %g%% noise. Layouts selected using the original 5%% results.', ...
        nl,100*Data.noise_fractions(ib)), ...
        'All layouts shown; boxes: Q1-Q3; whiskers: 1.5 IQR; dark line connects median centers.'},Options);
    Figures.(figure_fields{k})=fig;
    study_export_figure(fig,output_dir,file_names{k},Options);
end

% Recall requested by the user is kept as a separate, readable main figure.
fig=make_figure('Recall versus noise',[120 100 1150 760],Options);
ax=axes('Parent',fig,'Position',[.095 .23 .86 .63]);hold(ax,'on');
[mu,sd]=metric_summary(Summary,Data,'Recall');
draw_noise_curves(ax,100*Data.noise_fractions,100*mu,100*sd,colors,Options,true);
set(ax,'YLim',[0 105],'YTick',0:20:100);ylabel(ax,'Mean recall (%)');xlabel(ax,'Nominal noise level (% of clean-map RMS)');
legend(ax,labels,'Location','southwest','Box','off','FontName',Options.font_name,...
    'FontSize',Options.font_size);style_axes(ax,Options);
heading(fig,sprintf('Recall versus noise | same %d selected layouts at every noise level',nl),Options);
if Options.show_sd
    spread='Lines: layout means; shading: mean +/- 1 sample SD (not a confidence interval; display clipped to 0-100%).';
else,spread='Lines: layout means. The full per-layout values and sample SD are saved in the numerical tables.';end
note(fig,{spread,'Within each layout, the same standard Gaussian noise map is scaled across noise levels; noise-level results are paired.',...
    'The layouts were selected using results at 5% noise. The other levels do not form a new independent layout sample.'},Options);
Figures.recall_noise=fig;
study_export_figure(fig,output_dir,'Ensemble02_recall_vs_noise',Options);

% Supporting panels retain the trade-off that a recall-only curve may hide.
fig=make_figure('Detection trade-offs versus noise',[120 100 1450 610],Options);
positions=[.055 .27 .27 .56;.385 .27 .27 .56;.715 .27 .27 .56];
fields={'F1','FP','GOSPA_m'};scales=[100 1 1];
titles={'(a) F1 score','(b) Unmatched outputs','(c) Standard GOSPA'};
ylabels={'Mean F1 (%)','Mean false positives per map','Mean GOSPA (m)'};
for k=1:3
    ax=axes('Parent',fig,'Position',positions(k,:));hold(ax,'on');
    [mu,sd]=metric_summary(Summary,Data,fields{k});
    draw_noise_curves(ax,100*Data.noise_fractions,scales(k)*mu,scales(k)*sd,colors,Options,k==1);
    if k==1,set(ax,'YLim',[0 105],'YTick',0:20:100);else,ylim(ax,[0 max(1,max(mu(:)+sd(:)))*1.08]);end
    title(ax,titles{k},'FontWeight','bold');ylabel(ax,ylabels{k});xlabel(ax,'Nominal noise (%)');
    style_axes(ax,Options);
    if k==1,legend(ax,labels,'Location','southwest','Box','off',...
            'FontName',Options.font_name,'FontSize',Options.font_size-1);end
end
heading(fig,sprintf('Supporting metrics | same %d selected layouts | no pooling across noise levels',nl),Options);
if Options.show_sd,spread='Lines: per-noise means over layouts; shading: mean +/- 1 sample SD, not a confidence interval.';
else,spread='Lines: per-noise means over layouts; sample SD and individual results are saved in the numerical tables.';end
note(fig,{spread,'Fixed layouts were selected at 5% noise. The same noise pattern is scaled within each layout across the noise sweep.'},Options);
Figures.tradeoff_noise=fig;
study_export_figure(fig,output_dir,'Ensemble03_supporting_noise_metrics',Options);
end

function Options=defaults(Options)
d=struct('visible','on','save_figures',true,'save_editable_fig',false,...
    'font_name','Times New Roman','font_size',12,'dpi',300,'box_noise_fraction',.05,...
    'show_sd',true,'pair_lines',false,'show_titles',false,'colors',[.36 .44 .57;.84 .55 .19;.03 .49 .55]);
fields=fieldnames(d);for k=1:numel(fields),if ~isfield(Options,fields{k}),Options.(fields{k})=d.(fields{k});end,end
assert(isequal(size(Options.colors),[3 3])&&all(Options.colors(:)>=0)&&all(Options.colors(:)<=1));
end

function fig=make_figure(name,position,Options)
fig=figure('Name',name,'Color','w','Position',position,'Visible',Options.visible,'Renderer','painters');
end

function style_axes(ax,Options)
set(ax,'FontName',Options.font_name,'FontSize',Options.font_size,'FontWeight','bold','LineWidth',1.5,...
    'TickDir','in','Box','on','Layer','top','XGrid','off','YGrid','off',...
    'GridColor',[.82 .82 .82],'GridAlpha',.35);
end

function heading(fig,text,Options)
annotation(fig,'textbox',[.03 .92 .94 .06],'String',text,'EdgeColor','none',...
    'HorizontalAlignment','center','VerticalAlignment','middle','FontName',Options.font_name,...
    'FontSize',Options.font_size+3,'FontWeight','bold','Interpreter','none');
end

function note(fig,text,Options)
% Separate single-line annotations also work with Octave's gnuplot backend.
height=.12/numel(text);
for k=1:numel(text)
    y=.015+.12-k*height;
    annotation(fig,'textbox',[.04 y .92 height],'String',text{k},'EdgeColor','none',...
        'HorizontalAlignment','center','VerticalAlignment','middle','FontName',Options.font_name,...
        'FontSize',max(8,Options.font_size-2),'Interpreter','none');
end
end

function draw_paired_boxes(ax,values,colors,pair_lines)
n=size(values,1);offsets=linspace(-.075,.075,n).';
if pair_lines
    for i=1:n,plot(ax,(1:3)+offsets(i),values(i,:),'-','Color',[.80 .80 .80],'LineWidth',.6,'HandleVisibility','off');end
end
for j=1:3
    v=values(:,j);assert(all(isfinite(v)),'Boxplots require finite primary metrics for all layouts.');
    quartiles=linear_quantile(v,[.25 .5 .75]);q1=quartiles(1);md=quartiles(2);q3=quartiles(3);iqr=q3-q1;
    inside=v>=q1-1.5*iqr & v<=q3+1.5*iqr;
    low=min(v(inside));high=max(v(inside));w=.42;
    fillcolor=.84+.16*colors(j,:);
    patch(ax,j+[-w/2 w/2 w/2 -w/2],[q1 q1 q3 q3],fillcolor,...
        'EdgeColor',colors(j,:),'LineWidth',1.6,'HandleVisibility','off');
    plot(ax,[j j],[low q1],'-','Color',colors(j,:),'LineWidth',1.2,'HandleVisibility','off');
    plot(ax,[j j],[q3 high],'-','Color',colors(j,:),'LineWidth',1.2,'HandleVisibility','off');
    plot(ax,j+[-.12 .12],[low low],'-','Color',colors(j,:),'LineWidth',1.2,'HandleVisibility','off');
    plot(ax,j+[-.12 .12],[high high],'-','Color',colors(j,:),'LineWidth',1.2,'HandleVisibility','off');
    plot(ax,j+[-w/2 w/2],[md md],'-','Color',[.10 .10 .10],'LineWidth',1.8,'HandleVisibility','off');
    plot(ax,j+offsets,v,'o','Color',colors(j,:),'MarkerFaceColor',colors(j,:),...
        'MarkerSize',4,'LineStyle','none','HandleVisibility','off');
end
% Draw last so the line reaches each median center instead of vanishing
% behind the filled boxes. All per-layout scatter points remain present.
plot(ax,1:3,median(values,1),'-o','Color',[.2 .2 .2],'MarkerFaceColor',[.2 .2 .2],...
    'MarkerSize',4,'LineWidth',1.4,'HandleVisibility','off');
end

function q=linear_quantile(values,p)
values=sort(values(:));h=1+(numel(values)-1)*p;lo=floor(h);hi=ceil(h);
q=values(lo).'+(h-lo).*(values(hi).'-values(lo).');
end

function [mu,sd]=metric_summary(Summary,Data,field)
nn=numel(Data.noise_fractions);mu=NaN(nn,3);sd=mu;
for i=1:nn
    for j=1:3
        found=[];
        for k=1:numel(Summary)
            if abs(Summary(k).noise_fraction-Data.noise_fractions(i))<=1e-12&&strcmp(Summary(k).method,Data.method_names{j})
                found(end+1)=k; %#ok<AGROW>
            end
        end
        assert(numel(found)==1,'Summary lacks one unique method/noise combination.');q=Summary(found);
        assert(q.n_layouts==numel(Data.layout_ids),'Summary is not based on all selected layouts.');
        mu(i,j)=q.([field '_mean']);sd(i,j)=q.([field '_SD']);
    end
end
end

function draw_noise_curves(ax,x,mu,sd,colors,Options,percent)
x=x(:).';
if Options.show_sd
    for j=1:3
        lo=max(0,mu(:,j)-sd(:,j));hi=mu(:,j)+sd(:,j);
        if percent,hi=min(100,hi);end
        shade=.88+.12*colors(j,:);
        fill(ax,[x fliplr(x)],[lo(:).' fliplr(hi(:).')],shade,...
            'EdgeColor','none','HandleVisibility','off');
    end
end
markers={'o','s','^'};
for j=1:3
    plot(ax,x,mu(:,j),['-' markers{j}],'Color',colors(j,:),'MarkerFaceColor',colors(j,:),...
        'LineWidth',1.8,'MarkerSize',6);
end
set(ax,'XTick',x);
if numel(x)>1,pad=max(.25,.025*(max(x)-min(x)));xlim(ax,[min(x)-pad max(x)+pad]);
else,xlim(ax,[x-.5 x+.5]);end
end
