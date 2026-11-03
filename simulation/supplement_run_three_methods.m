function R=supplement_run_three_methods(Scene,SourceInfo,Meta,S)
% Fresh three-method pipeline using the archived numerical core unchanged.
% The unused historical pairwise game is omitted. Its local-cost initial
% actions are exactly the same as State.initial_actions. No cache is read.
started=tic;Data=struct('X',Scene.X,'Y',Scene.Y,'Bz',Scene.Bobserved);
t=tic;Pool=pg_pool(Data.Bz,Data.X,Data.Y,S.Candidate);candidate_seconds=toc(t);
fprintf('Candidate pool: %d branches, %d nodes, %.1f s\n', ...
    Pool.branch_count,Pool.node_count,candidate_seconds);
[PriminFit,PriminInfo]=mw_primin_baseline(Data,Pool,S.Primin);
fprintf('Primin baseline: %d estimates\n',size(PriminFit.position,1));
[Windows,LocalInfo]=mw_local_invert(Data,Pool,S.Local);
t=tic;State=mw_physical_setup(Data,Windows,S.Fusion,S.Physical);setup_seconds=toc(t);
initial=State.initial_actions;
[ind,ind_pi,ind_fi]=mw_physical_score(initial,State);
[actions,GI,phy,phy_pi,phy_fi]=mw_physical_game_fusedcache(State,S.Game);
assert(isequal(GI.initial_actions,initial));
if ~GI.converged
    warning('WindowGame:SweepLimit', ...
        'Sweep limit reached before unilateral stability; this is recorded in GameInfo.');
end
Fits={PriminFit,ind,phy};Names={'primin_direct','independent_refit','physical_game'};
Actions={[],initial,actions};PIs={struct(),ind_pi,phy_pi};
FIs={struct(),ind_fi,phy_fi};GIs={struct(),struct(),GI};
Methods=cell(3,1);small=cell(3,1);
for j=1:3
    [M,D,E]=multi_evaluate_output(Fits{j},Scene,PIs{j},j>1,S.Evaluation.position_tolerance_m);
    Methods{j}=struct('name',Names{j},'Fit',Fits{j},'actions',Actions{j}, ...
        'Metrics',M,'Details',D,'Evaluation',E,'PhysicalInfo',PIs{j}, ...
        'GameInfo',GIs{j},'FusionInfo',FIs{j});
    small{j}=struct('name',Names{j},'Metrics',M,'Evaluation',E, ...
        'positions',Fits{j}.position,'moments',Fits{j}.moment);
    fprintf('%s: TP %d FP %d FN %d F1 %.4f GOSPA %.6g m\n', ...
        Names{j},M.TP,M.FP,M.FN,M.F1,M.GOSPA);
end
Summary=Meta;Summary.Methods=small;Summary.Truth=Scene.Truth;
Summary.SourceInfo=SourceInfo;Summary.sigma_hat_nT=State.sigma_hat_nT;
Summary.actual_injected_sigma_nT=Scene.Sigma_nT;Summary.tau=State.tau;
Summary.elapsed_seconds=toc(started);
Summary.Diagnostics=struct('candidate_seconds',candidate_seconds, ...
    'primin_seconds',PriminInfo.elapsed_seconds, ...
    'local_seconds',LocalInfo.elapsed_seconds,'physical_setup_seconds',setup_seconds, ...
    'physical_selection_seconds',GI.selection_seconds, ...
    'accepted_nonconverged_hypotheses',LocalInfo.accepted_nonconverged_hypotheses);
R=struct('Meta',Meta,'Scene',Scene,'SourceInfo',SourceInfo,'Methods',{Methods}, ...
    'Summary',Summary,'elapsed_seconds',Summary.elapsed_seconds);
R.Computation=struct('physical_engine','mw_physical_game_fusedcache', ...
    'mode','fresh recomputation','software',version,'computer',computer, ...
    'note','No cached estimates, local fits or candidates were reused.');
R=wg_compact_result(R);
end
