function res=CropImg(img,nxnew,nynew);
% take the central part of a 2D image
nx=size(img,1);
ny=size(img,2);
N1=(nx-nxnew)/2+1:nx-(nx-nxnew)/2;
N2=(ny-nynew)/2+1:ny-(ny-nynew)/2;
res=img(N1,N2,:,:,:,:,:);

return


