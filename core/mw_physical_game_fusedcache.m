function [actions,Info,Fit,PhysicalInfo,FusionInfo]=mw_physical_game_fusedcache(State,Game)
%MW_PHYSICAL_GAME Finite common-interest game for the actual fused output.
% Each action is a fixed local source-set hypothesis.  u_i=-J makes -J an
% exact potential.  This is also discrete coordinate descent on the same J.
% Accepted moves strictly decrease J; all final unilateral alternatives are
% explicitly checked.  No true source information enters the update.
if nargin<2,Game=struct();end
if ~isfield(Game,'max_sweeps'),Game.max_sweeps=50;end
if ~isfield(Game,'improvement_tolerance'),Game.improvement_tolerance=1e-10;end
if ~isfield(Game,'order'),Game.order='forward';end
if ~isfield(Game,'verbose'),Game.verbose=false;end
assert(Game.max_sweeps>=1&&Game.max_sweeps==round(Game.max_sweeps));
assert(Game.improvement_tolerance>=0&&isfinite(Game.improvement_tolerance));
assert(any(strcmp(Game.order,{'forward','reverse'})));
started=tic;n=numel(State.Windows);actions=State.initial_actions;
cache=containers.Map('KeyType','char','ValueType','double');
fused_cache=containers.Map('KeyType','char','ValueType','double');
fused_evaluations=0;fused_cache_hits=0;cache_validation_checks=0;
evaluations=0;cache_hits=0;
current=score(actions);initial_energy=current;energy_trace=current;
changes_per_sweep=zeros(0,1);
template=struct('sweep',0,'window',0,'old_action',0,'new_action',0,...
 'old_energy',NaN,'new_energy',NaN,'energy_change',NaN);
accepted=repmat(template,0,1);
if strcmp(Game.order,'reverse'),visit=n:-1:1;else,visit=1:n;end
for sweep=1:Game.max_sweeps
 changed=0;
 for i=visit
  old_action=actions(i);best=current;best_action=old_action;
  for a=State.valid_actions{i}
   if a==old_action,continue;end
   trial=actions;trial(i)=a;value=score(trial);
   if value<best,best=value;best_action=a;end
  end
  threshold=Game.improvement_tolerance*max(1,abs(current));
  if best<current-threshold
   rec=template;rec.sweep=sweep;rec.window=i;rec.old_action=old_action;
   rec.new_action=best_action;rec.old_energy=current;rec.new_energy=best;
   rec.energy_change=best-current;accepted(end+1,1)=rec; %#ok<AGROW>
   actions(i)=best_action;current=best;changed=changed+1;
  end
 end
 energy_trace(end+1,1)=current;changes_per_sweep(end+1,1)=changed; %#ok<AGROW>
 if Game.verbose,fprintf('  Physical sweep %d: J=%.9g, %d changes\n',sweep,current,changed);end
 if changed==0,break;end
end
max_gain=0;terminal_best_actions=actions;
for i=1:n
 best=current;
 for a=State.valid_actions{i}
  trial=actions;trial(i)=a;value=score(trial);
  if value<best,best=value;terminal_best_actions(i)=a;end
 end
 max_gain=max(max_gain,current-best);
end
[Fit,PhysicalInfo,FusionInfo]=mw_physical_score_fast(actions,State);
assert(abs(PhysicalInfo.J-current)<=1e-10*max(1,abs(current)),...
 'The same action profile received a different deterministic score.');
assert(all(diff(energy_trace)<=0)&&current<=initial_energy);
terminal_tolerance=Game.improvement_tolerance*max(1,abs(current));
Info=struct('initial_actions',State.initial_actions,'initial_energy',initial_energy,...
 'final_energy',current,'energy_trace',energy_trace,'potential_trace',-energy_trace,...
 'changes_per_sweep',changes_per_sweep,'changes',numel(accepted),...
 'changed_windows',sum(actions~=State.initial_actions),...
 'sweeps',numel(changes_per_sweep),'max_gain',max_gain,...
 'terminal_tolerance',terminal_tolerance,'converged',max_gain<=terminal_tolerance,...
 'terminal_best_actions',terminal_best_actions,'accepted_moves',accepted,...
 'score_evaluations',evaluations,'cache_hits',cache_hits,'order',Game.order,...
 'selection_seconds',toc(started),'tau',State.tau,...
 'objective','global fixed-position refit SSE/(2*sigma_hat^2) + tau*K');
Info.evaluations=Info.score_evaluations;
Info.fused_score_evaluations=fused_evaluations;
Info.fused_cache_hits=fused_cache_hits;
Info.fused_cache_validation_checks=cache_validation_checks;

 function value=score(a)
  key=sprintf('%d,',a);
  if isKey(cache,key)
   value=cache(key);cache_hits=cache_hits+1;
  else
   % Numerical optimization only: the ordered fused representatives determine
   % the exact same global SVD design. No score, candidate, or fusion changes.
   f=mw_fuse_fast(State.Windows,a,State.Data,State.Fusion);
   ids=zeros(0,1);
   for wi=1:numel(a),ids=[ids;State.report_ids{wi}{a(wi)}(:)];end
   if isempty(f.position),gids=zeros(0,1);else
    reps=reshape([f.clusters.representative_report_id],[],1);gids=ids(reps);
   end
   fusedkey=['f:',sprintf('%d,',gids)];
   if isKey(fused_cache,fusedkey)
    value=fused_cache(fusedkey);fused_cache_hits=fused_cache_hits+1;
    if cache_validation_checks<10
     [~,p]=mw_physical_score(a,State);
     assert(abs(p.J-value)<=1e-12*max(1,abs(value)));
     cache_validation_checks=cache_validation_checks+1;
    end
   else
    [~,p]=mw_physical_score_fast(a,State);value=p.J;
    fused_cache(fusedkey)=value;fused_evaluations=fused_evaluations+1;
   end
   assert(isfinite(value));cache(key)=value;evaluations=evaluations+1;
  end
 end
end
