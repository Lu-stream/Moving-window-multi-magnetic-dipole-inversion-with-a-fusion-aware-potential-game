function study_export_figure(f,out,name,O)
%STUDY_EXPORT_FIGURE Export full-figure PNG/PDF, optionally editable figures.
if nargin<4,O=struct();end
if ~isfield(O,'save_figures'),O.save_figures=true;end
if ~isfield(O,'save_editable_fig'),O.save_editable_fig=false;end
if ~isfield(O,'dpi'),O.dpi=300;end
if ~exist(out,'dir'),mkdir(out);end
if ~O.save_figures,return;end
set(f,'PaperPositionMode','auto','InvertHardcopy','off');
pixel_position=getpixelposition(f);page_size=pixel_position(3:4)/100;
set(f,'PaperUnits','inches','PaperPosition',[0 0 page_size],'PaperSize',page_size);
base=fullfile(out,name);
if exist('OCTAVE_VERSION','builtin')
    % Stage new bytes so renderer cleanup cannot truncate a finished export.
    stage_base=tempname(tempdir);
    png_stage=[stage_base '.png'];pdf_stage=[stage_base '.pdf'];
    if exist(png_stage,'file'),delete(png_stage);end
    if exist(pdf_stage,'file'),delete(pdf_stage);end
    print(f,png_stage,'-dpng',sprintf('-r%d',O.dpi));
    wait_png_complete(png_stage);
    print(f,pdf_stage,'-dpdf');
    wait_pdf_complete(pdf_stage);
    movefile(png_stage,[base '.png'],'f');
    movefile(pdf_stage,[base '.pdf'],'f');
    if O.save_editable_fig,hgsave(f,[base '.ofig']);end
else
    % print supports the complete figure (including annotations) in R2021b.
    pixel_position=getpixelposition(f);page_size=pixel_position(3:4)/100;
    set(f,'PaperUnits','inches','PaperPosition',[0 0 page_size],'PaperSize',page_size);
    print(f,[base '.png'],'-dpng',sprintf('-r%d',O.dpi));
    print(f,[base '.pdf'],'-dpdf','-painters');
    if O.save_editable_fig,savefig(f,[base '.fig']);end
end
end

function wait_png_complete(path)
% GNUPlot may finish its output pipe asynchronously. Do not close a figure
% until its PNG has an IEND chunk; MATLAB native print does not use this path.
for attempt=1:50
    fid=fopen(path,'rb');complete=false;
    if fid>=0
        if fseek(fid,-12,'eof')==0
            tail=fread(fid,12,'uint8').';
            complete=numel(tail)==12&&isequal(tail(5:8),double('IEND'));
        end
        fclose(fid);
    end
    if complete,return;end
    pause(.1);
end
error('PNG export did not complete: %s',path);
end

function wait_pdf_complete(path)
for attempt=1:50
    fid=fopen(path,'rb');complete=false;
    if fid>=0
        fseek(fid,0,'eof');n=ftell(fid);fseek(fid,max(0,n-512),'bof');
        tail=char(fread(fid,Inf,'uint8').');fclose(fid);
        complete=~isempty(strfind(tail,'%%EOF')); %#ok<STREMP>
    end
    if complete,return;end
    pause(.1);
end
error('PDF export did not complete: %s',path);
end
