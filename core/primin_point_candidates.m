function out=primin_point_candidates(Bz,G,r_obs)
%PRIMIN_POINT_CANDIDATES Analytic dipole candidates adapted from text_inverse.m.
% Source: user's elsarticle.rar / CG投稿文献/代码可用性/text_inverse.m.
% Theory: Primin & Nedayvoda (2023), doi:10.1007/s10559-023-00618-7.
% Bz nT; G nT/m; observation Cartesian coordinates use z POSITIVE UP.
% Returned position = [x,y,depth below observation plane], depth positive.
% The original lambda/n/normalized-moment/r formulas are retained.
% We permute INDICES, not eigenvalue values matched with an absolute tolerance.
    out=struct('position',zeros(0,3),'moment',zeros(0,3), ...
        'self_error',zeros(0,1),'lambda',nan(1,3));
    if ~all(isfinite(G(:))) || ~isfinite(Bz) || norm(G,'fro')==0, return; end
    G=(G+G.')/2;
    [V,D]=eig(G); lambda=diag(D);
    orders=perms(1:3); valid=false(size(orders,1),1); scores=-inf(size(valid));
    for k=1:size(orders,1)
        a=lambda(orders(k,1)); b=lambda(orders(k,2));
        if abs(b-a)<=1e-12*max(abs(lambda)), continue; end
        s2=(2*b+a)/(b-a);
        valid(k)=(s2>=-1e-10 && s2<=1+1e-10);
        if valid(k), scores(k)=a+2*b; end % original permutation ranking
    end
    if ~any(valid), return; end
    [~,best]=max(scores); ix=orders(best,:);
    lam=lambda(ix); E=V(:,ix); out.lambda=lam(:).';
    s2=max(0,min(1,(2*lam(2)+lam(1))/(lam(2)-lam(1))));
    sin_values=sqrt(s2)*[1 -1 1 -1];
    cos_values=sqrt(1-s2)*[1 1 -1 -1];
    for k=1:4
        n1=cos_values(k); n2=sin_values(k);
        n=E*[n1;n2;0];
        normalized=E*[lam(1)*n1*(5*n2^2-1)/(1+n1^2); ...
            lam(1)*n2*(1-5*n1^2)/(1+n1^2);0];
        denominator=3*dot(normalized,n)*n(3)-normalized(3);
        if abs(denominator)<=1e-12*norm(normalized), continue; end
        distance=3*Bz/denominator;
        % Distance must be positive and the source below the observation plane.
        if ~isfinite(distance) || distance<=0 || n(3)<=0, continue; end
        r=distance*n;
        pos=[r_obs(1)-r(1),r_obs(2)-r(2),r(3)];
        moment=(distance^4/(3*100))*normalized;
        Bcheck=100*(3*dot(moment,r)*r/norm(r)^5-moment/norm(r)^3);
        out.position(end+1,:)=pos;
        out.moment(end+1,:)=moment(:).';
        out.self_error(end+1,1)=abs(Bcheck(3)-Bz);
    end
end
