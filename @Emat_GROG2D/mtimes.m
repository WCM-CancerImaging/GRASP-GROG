function res = mtimes(a,b)

if a.adjoint
    % to image space: n-coil -> 1-coil
    b=b.*a.mask.*a.W1;
    res=sum(ifft2c_mri(b).*conj(a.b1),4)./sum(abs((a.b1)).^2,4);
else
    % to k-space: 1-coil -> n-coil
    res=fft2c_mri(repmat(b,[1,1,1,size(a.b1,4)]).*a.b1).*a.mask;
    res=res.*a.W1;
end





    
