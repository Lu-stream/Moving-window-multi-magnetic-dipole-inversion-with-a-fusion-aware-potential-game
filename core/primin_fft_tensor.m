function tensor=primin_fft_tensor(Bz,dx,dy,pad_factor)
%PRIMIN_FFT_TENSOR FFT reconstruction above sources, z positive upward.
% Rows are y; columns are x; no output transposition.
% Corrects Fourier multipliers/axis indexing in the uploaded Fuliye_jisuan.m.
% Uses zero padding; finite survey truncation still introduces errors.
    [ny,nx]=size(Bz); py=ceil(pad_factor*ny); px=ceil(pad_factor*nx);
    assert(py>=ny && px>=nx && dx>0 && dy>0);
    sy=floor((py-ny)/2); sx=floor((px-nx)/2);
    iy=sy+(1:ny); ix=sx+(1:nx);
    padded=zeros(py,px); padded(iy,ix)=Bz;
    kx=2*pi/(px*dx)*[0:ceil(px/2)-1,-floor(px/2):-1];
    ky=2*pi/(py*dy)*[0:ceil(py/2)-1,-floor(py/2):-1];
    [KX,KY]=meshgrid(kx,ky); K=hypot(KX,KY);
    safe=K; safe(K==0)=Inf;
    % Odd derivatives at the self-conjugate Nyquist bins must be zero.
    kxd=kx; kyd=ky;
    if mod(px,2)==0, kxd(px/2+1)=0; end
    if mod(py,2)==0, kyd(py/2+1)=0; end
    [DX,DY]=meshgrid(kxd,kyd);
    F=fft2(padded);
    multipliers={KX.^2./safe,DX.*DY./safe,1i*DX, ...
        KY.^2./safe,1i*DY,-K};
    labels={'xx','xy','xz','yy','yz','zz'};
    tensor=struct();
    for k=1:numel(labels)
        value=real(ifft2(F.*multipliers{k}));
        tensor.(labels{k})=value(iy,ix);
    end
end
