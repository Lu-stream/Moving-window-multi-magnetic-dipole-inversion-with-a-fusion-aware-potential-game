function C=area1_prepare_candidates(input_dir,output_dir,mode)
% Explicit candidate modes: missing inputs raise an error, never a fallback.
if strcmp(mode,'fixed_candidates')
 file=fullfile(input_dir,'area1_fixed_candidates.mat');
 assert(exist(file,'file')==2,'Missing fixed candidate file: %s',file);
 C=load(file);
 [Candidate,Primin,Local,Fusion]=area1_candidate_settings();
 check_settings(C.Candidate,Candidate,'Candidate');
 check_settings(C.Primin,Primin,'Primin');
 check_settings(C.Local,Local,'Local');
 check_settings(C.Fusion,Fusion,'Fusion');
 return;
end
assert(strcmp(mode,'rebuild_candidates'),'Unknown candidate mode: %s',mode);
raw_file=fullfile(input_dir,'Area1-90-fig2MMT.txt');
assert(exist(raw_file,'file')==2,'Missing raw magnetic map: %s',raw_file);
Bz=dlmread(raw_file)*1e9;
assert(isequal(size(Bz),[201 351]),'Expected the registered 201-by-351 Area 1 map.');
[X,Y]=meshgrid(-175:175,-100:100);
Data=struct('Bz',Bz,'X',X,'Y',Y);
[Candidate,Primin,Local,Fusion]=area1_candidate_settings();
% Source search derives from scan extent and the declared 0--74 um depth
% below the sample, plus the known 2 um sensor-to-surface distance.
Candidate.window_widths_m=1.001*max([max(X(:))-min(X(:)),max(Y(:))-min(Y(:))]);
fprintf('Rebuilding the analytical pool from the raw map...\n');
t=tic; Pool=pg_pool(Data.Bz,Data.X,Data.Y,Candidate);pool_seconds=toc(t);
save(fullfile(output_dir,'checkpoint_pool_rebuilt.mat'),'Data','Pool','Candidate','pool_seconds','-v7');
fprintf('Rebuilding Primin baseline (the scan may take several minutes)...\n');
[PriminFit,PriminInfo,PriminRaw]=mw_primin_baseline(Data,Pool,Primin);
save(fullfile(output_dir,'checkpoint_primin_rebuilt.mat'),'PriminFit','PriminInfo','PriminRaw','-v7');
fprintf('Rebuilding local nonlinear fits (potentially a long computation)...\n');
[Windows,LocalInfo]=mw_local_invert(Data,Pool,Local);
save(fullfile(output_dir,'area1_candidates_rebuilt.mat'),'Data','Candidate',...
 'Windows','LocalInfo','Local','PriminFit','PriminInfo','Fusion','Primin','-v7');
C=struct('Data',Data,'Candidate',Candidate,'Windows',Windows,...
 'LocalInfo',LocalInfo,'Local',Local,'PriminFit',PriminFit,...
 'PriminInfo',PriminInfo,'Fusion',Fusion,'Primin',Primin);
end

function check_settings(cached,current,name)
fields=fieldnames(current);
for k=1:numel(fields)
 f=fields{k};assert(isfield(cached,f),'Fixed cache has no setting %s.%s.',name,f);
 a=cached.(f);b=current.(f);
 if isnumeric(b)||islogical(b)
  same=isequal(size(a),size(b))&&all(abs(double(a(:))-double(b(:)))<=1e-12*max(1,max(abs(double(b(:))))));
 else,same=isequal(a,b);end
 assert(same,'Fixed candidate cache differs at %s.%s. Choose rebuild_candidates explicitly.',name,f);
end
end
