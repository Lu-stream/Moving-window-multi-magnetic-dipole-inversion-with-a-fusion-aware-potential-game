function [Metrics,Details]=pilot_evaluate(Fit,Truth,c)
%PILOT_EVALUATE Exact gated matching and GOSPA without external toolboxes.
% Truth.position and Fit.position are N-by-3 arrays in the same length unit.
% TP matching first maximizes the number of pairs with distance <= c, then
% minimizes their total distance. GOSPA uses p=1, alpha=2 and cutoff c.
% A subset dynamic program is exact for up to 16 true sources and arbitrary
% output count. The preliminary experiment has only two true sources.
% mean_matched_position_error concerns TP pairs only (NaN if none).
assert(isscalar(c)&&isfinite(c)&&c>0,'c must be finite and positive.');
P=positions(Truth);Q=positions(Fit);
K=size(P,1);Kh=size(Q,1);
assert(K<=16,'pilot_evaluate supports at most 16 true sources.');
dist=zeros(K,Kh);
for k=1:K
    dist(k,:)=sqrt(sum(bsxfun(@minus,Q,P(k,:)).^2,2)).';
end
states=2^K;
card=zeros(1,states);
for s=0:states-1
    card(s+1)=sum(bitget(s,1:K));
end
% dp(s) is the smallest total distance matching exactly the truth subset s.
% At each output point, either skip it or match it to one unused true source.
dp=Inf(1,states);dp(1)=0;
parent=zeros(Kh,states,'uint32');
chosen=zeros(Kh,states,'uint8');
for j=1:Kh
    next=dp;
    parent(j,:)=uint32(0:states-1);
    for idx=find(isfinite(dp))
        mask=idx-1;
        for k=1:K
            if ~bitget(mask,k)&&dist(k,j)<=c
                newmask=bitset(mask,k,1);newidx=newmask+1;
                value=dp(idx)+dist(k,j);
                if value<next(newidx)
                    next(newidx)=value;
                    parent(j,newidx)=uint32(mask);
                    chosen(j,newidx)=uint8(k);
                end
            end
        end
    end
    dp=next;
end
finite=find(isfinite(dp));
bestcard=max(card(finite));
eligible=finite(card(finite)==bestcard);
[~,ii]=min(dp(eligible));tpidx=eligible(ii);
[truth_ids,estimate_ids]=trace_assignment(tpidx-1,parent,chosen);
errors=zeros(numel(truth_ids),1);
for n=1:numel(truth_ids)
    errors(n)=dist(truth_ids(n),estimate_ids(n));
end
% For p=1, alpha=2, every associated pair replaces two c/2 penalties:
% G = c/2*(K+Kh) + min_partial_matching sum(distance-c).
% For a fixed truth subset, minimizing distance gives the same dp above.
[saving,gi]=min(dp-c*card);
gospa=max(0,c/2*(K+Kh)+saving);
[gtruth,gestimate]=trace_assignment(gi-1,parent,chosen);
TP=numel(truth_ids);FP=Kh-TP;FN=K-TP;
if K>0,recall=TP/K;else,recall=NaN;end
if Kh>0,precision=TP/Kh;else,precision=NaN;end
if K+Kh>0,f1=2*TP/(K+Kh);else,f1=1;end
if isempty(errors),meanerror=NaN;else,meanerror=mean(errors);end
Metrics=struct('TP',TP,'FP',FP,'FN',FN,'Recall',recall,...
    'Precision',precision,'F1',f1,'GOSPA',gospa,...
    'output_count',Kh,'mean_matched_position_error',meanerror);
Details=struct('matched_truth_indices',truth_ids,...
    'matched_estimate_indices',estimate_ids,'matched_errors',errors,...
    'unmatched_truth_indices',setdiff((1:K).',truth_ids),...
    'unmatched_estimate_indices',setdiff((1:Kh).',estimate_ids),...
    'distance_matrix',dist,'gospa_truth_indices',gtruth,...
    'gospa_estimate_indices',gestimate,'cutoff',c,...
    'empty_output_gospa',c*K/2);
end

function P=positions(S)
assert(isstruct(S)&&isfield(S,'position'),'Input requires .position.');
P=S.position;
if isempty(P),P=zeros(0,3);end
assert(isnumeric(P)&&ismatrix(P)&&size(P,2)==3,...
    'Positions must be N-by-3.');
assert(all(isfinite(P(:))),'Positions must be finite.');
end

function [truth_ids,estimate_ids]=trace_assignment(mask,parent,chosen)
truth_ids=zeros(0,1);estimate_ids=zeros(0,1);
for j=size(parent,1):-1:1
    k=double(chosen(j,mask+1));
    prev=double(parent(j,mask+1));
    if k>0
        truth_ids(end+1,1)=k; %#ok<AGROW>
        estimate_ids(end+1,1)=j; %#ok<AGROW>
    end
    mask=prev;
end
[truth_ids,order]=sort(truth_ids);
estimate_ids=estimate_ids(order);
end
