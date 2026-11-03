function yes=wg_same_scene(A,B)
% Cross-runtime last-bit differences must not invalidate a frozen-input replay.
% This tolerance is only for machine roundoff, never for scientific changes.
yes=isequal(size(A.X),size(B.X))&&isequal(size(A.Y),size(B.Y))&&...
 isequal(size(A.Bobserved),size(B.Bobserved))&&...
 isequal(A.Truth.position,B.Truth.position)&&...
 A.NoiseFraction==B.NoiseFraction&&A.NoiseSeed==B.NoiseSeed;
if ~yes,return;end
fields={'X','Y','Bclean','Bobserved','Noise'};
for q=1:numel(fields)
 a=A.(fields{q});b=B.(fields{q});
 if any(~isfinite(a(:)))||any(~isfinite(b(:)))||...
  max(abs(a(:)-b(:)))>1e-12*max(1,max(abs(b(:))))
  yes=false;return;
 end
end
a=A.Truth.moment;b=B.Truth.moment;
yes=isequal(size(a),size(b))&&...
 max(abs(a(:)-b(:)))<=1e-12*max(1,max(abs(b(:))));
end
