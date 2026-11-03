function Results=area1_run_comparison(input_dir,output_dir,candidate_mode,Options)
% Run three measured-data methods. No CT coordinates enter inference.
% candidate_mode is exactly 'fixed_candidates' or 'rebuild_candidates'.
if ~exist(output_dir,'dir'),mkdir(output_dir);end
tau_multiplier=Options.tau_multiplier;
svd_relative_tolerance=Options.svd_relative_tolerance;
max_sweeps=Options.max_sweeps;
evaluation_tolerance_um=Options.evaluation_tolerance_um;
tolerance_grid_um=Options.tolerance_grid_um;
started=tic;
C=area1_prepare_candidates(input_dir,output_dir,candidate_mode);
Data=C.Data; Candidate=C.Candidate; Windows=C.Windows; LocalInfo=C.LocalInfo;
Local=C.Local; PriminFit=C.PriminFit; PriminInfo=C.PriminInfo; Fusion=C.Fusion; Primin=C.Primin;
% Confirm cached measured field is the actual raw measurement, converted T->nT.
raw=dlmread(fullfile(input_dir,'Area1-90-fig2MMT.txt'));
assert(isequal(size(raw),size(Data.Bz)));
assert(max(abs(raw(:)*1e9-Data.Bz(:)))<1e-8);
assert(Local.max_iterations==1000&&Local.max_function_evaluations==4000);
Physical=struct('tau_multiplier',tau_multiplier,...
 'svd_relative_tolerance',svd_relative_tolerance,...
 'version','fusion-aware-physical-v1');
Game=struct('max_sweeps',max_sweeps,'improvement_tolerance',1e-10,...
 'order','forward','verbose',true);
fprintf('Area1 measured map: %dx%d; %d fixed windows; %d accepted hypotheses.\n',...
 size(Data.Bz,1),size(Data.Bz,2),numel(Windows),LocalInfo.accepted_hypotheses);
fprintf('Local optimization nonconverged accepted hypotheses: %d\n',...
 LocalInfo.accepted_nonconverged_hypotheses);
State=mw_physical_setup(Data,Windows,Fusion,Physical);
fprintf('Setup: %d reports; sigma_hat=%.9g nT; tau=%.9g\n',...
 size(State.report_positions,1),State.sigma_hat_nT,State.tau);
initial_actions=State.initial_actions;
[IndependentFit,IndependentInfo,IndependentFusion]=mw_physical_score(initial_actions,State);
fprintf('Independent + fusion: K=%d RMSE=%.9g nT J=%.9g\n',...
 IndependentInfo.K,IndependentInfo.RMSE,IndependentInfo.J);
save(fullfile(output_dir,'checkpoint_independent.mat'),'IndependentFit',...
 'IndependentInfo','IndependentFusion','initial_actions','Physical','Game','-v7');
[game_actions,GameInfo,GameFit,GamePhysical,GameFusion]=mw_physical_game_fusedcache(State,Game);
% Independently run original (non-fast) fusion/score on final profile.
[CheckFit,CheckInfo]=mw_physical_score(game_actions,State);
assert(isequal(CheckFit.position,GameFit.position));
assert(abs(CheckInfo.J-GamePhysical.J)<=1e-10*max(1,abs(CheckInfo.J)));
if ~GameInfo.converged
 warning('area1:SweepLimit','Sweep limit reached; this output is not certified as unilaterally stable.');
end
assert(all(diff(GameInfo.energy_trace)<=0));
sigma_hat_nT=State.sigma_hat_nT;tau=State.tau;
clear State CheckFit CheckInfo;

% Primin direct reconstruction: same convention as latest synthetic baseline.
primin_prediction=zeros(size(Data.Bz));
for k=1:size(PriminFit.position,1)
 primin_prediction=primin_prediction+reshape(...
  dipole_bz_kernel(Data.X,Data.Y,PriminFit.position(k,:))*PriminFit.moment(k,:)',size(Data.Bz));
end
primin_residual=Data.Bz-primin_prediction;
PriminPhysical=struct('prediction',primin_prediction,'residual',primin_residual,...
 'RMSE',sqrt(mean(primin_residual(:).^2)),...
 'relative_residual',norm(primin_residual(:))/norm(Data.Bz(:)),...
 'J',NaN,'K',size(PriminFit.position,1));

%% Evaluation begins here: registered CT geometry has not entered inference.
xy_offset_um=[175 100];sensor_height_um=2;
Reference=area1_reference_from_cuboids(input_dir,output_dir,xy_offset_um,sensor_height_um);
Fits={PriminFit,IndependentFit,GameFit};
Names={'primin_direct','independent_refit','physical_game'};
Infos={PriminPhysical,IndependentInfo,GamePhysical};
Methods=cell(3,1);Summary=zeros(3,11);ToleranceResults=cell(3,numel(tolerance_grid_um));
for j=1:3
 F=Fits{j};[M,D]=pilot_evaluate(F,Reference,evaluation_tolerance_um);
 F.position_registered_um=bsxfun(@plus,F.position,[xy_offset_um 0]);
 F.moment_Am2=F.moment*1e-18;
 Methods{j}=struct('name',Names{j},'Fit',F,'Metrics',M,'Details',D,'PhysicalInfo',Infos{j});
 Summary(j,:)=[M.output_count M.TP M.FP M.FN M.Recall M.Precision M.F1 M.GOSPA,...
  M.mean_matched_position_error Infos{j}.RMSE Infos{j}.relative_residual];
 for t=1:numel(tolerance_grid_um)
  [TM,TD]=pilot_evaluate(F,Reference,tolerance_grid_um(t));
  ToleranceResults{j,t}=struct('Metrics',TM,'Details',TD);
 end
 fprintf('%s: K=%d CTmatches=%d unmatched=%d missed=%d recall=%.3f GOSPA=%.6g um RMSE=%.6g nT\n',...
  Names{j},M.output_count,M.TP,M.FP,M.FN,M.Recall,M.GOSPA,Infos{j}.RMSE);
end
Metadata=struct('data','Area1-90-fig2MMT.txt','source_doi','10.1594/PANGAEA.886724',...
 'measurement','SQUID/SSM random magnetic state','internal_length_unit','um',...
 'internal_field_unit','nT','internal_moment_to_Am2',1e-18,...
 'xy_offset_um',xy_offset_um,'sensor_height_um',sensor_height_um,...
 'evaluation_tolerance_um',evaluation_tolerance_um,'tolerance_grid_um',tolerance_grid_um,...
 'sigma_hat_nT',sigma_hat_nT,'tau',tau,'cached_local_candidates',strcmp(candidate_mode,'fixed_candidates'),...
 'candidate_mode',candidate_mode,'game_rerun',true,'engine','mw_physical_game_fusedcache','runtime_version',version,...
 'elapsed_seconds',toc(started),'run_timestamp_local',datestr(now,30));
Results=struct('Data',Data,'Reference',Reference,'Methods',{Methods},'Summary',Summary,...
 'GameInfo',GameInfo,'game_actions',game_actions,'initial_actions',initial_actions,...
 'LocalInfo',LocalInfo,'Local',Local,'Fusion',Fusion,'Physical',Physical,...
 'ToleranceResults',{ToleranceResults},'Metadata',Metadata);
save(fullfile(output_dir,'area1_latest_result.mat'),'Results','-v7');
fid=fopen(fullfile(output_dir,'summary.csv'),'w');
fprintf(fid,'method,outputs,CT_matches,unmatched_outputs,unmatched_CT,recall,precision,F1,GOSPA_um,mean_matched_error_um,RMSE_nT,relative_residual\n');
for j=1:3
 fprintf(fid,'%s',Names{j});fprintf(fid,',%.12g',Summary(j,:));fprintf(fid,'\n');
end
fclose(fid);
fid=fopen(fullfile(output_dir,'game_diagnostics.json'),'w');
fprintf(fid,'%s',jsonencode(struct('Metadata',Metadata,'GameInfo',GameInfo)));fclose(fid);
area1_write_supplementary(Results,output_dir);
fprintf('Completed in %.2f s. Saved %s\n',Metadata.elapsed_seconds,output_dir);

end
