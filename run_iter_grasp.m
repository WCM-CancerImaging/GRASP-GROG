function [grasp_recon, nufft_recon, b1] = run_iter_grasp(kdata,traj,dcf,spokes_per_frame,option)
cd 'F:\Research\GRASP_Breast\Simulation'

if strcmpi(option.method,'nufft')
    addpath('nufft_toolbox/');
    nx = size(kdata,1);
    ntviews = size(kdata,2);
%     disp('dcf- ramp-lak')
%     dcf = repmat(abs(linspace(-1,1,nx))',[1,ntviews]);

    nt = 1;
    nline = floor(ntviews/nt);
    for ii = 1:nt
        kdatau(:,:,:,ii) = kdata(:,(ii-1)*nline+1:ii*nline,:);
        dcfu(:,:,ii) = dcf(:,(ii-1)*nline+1:ii*nline);
        traju(:,:,ii) = traj(:,(ii-1)*nline+1:ii*nline);
    end
    param.E = MCNUFFT(traju,dcfu,ones(nx/2,nx/2));
%     param.y=kdatau.*repmat(reshape(dcfu,[size(dcfu,1),size(dcfu),1,1]),[1,1,1,size(kdatau,4)]);
    param.y=double(kdatau.*repmat(sqrt(dcfu),[1,1,size(kdatau,3),1]));
    clear b1 smap ref
    for c = 1:size(kdatau,3)
        ref(:,:,:,c) = param.E'*param.y(:,:,c,:);
    end
    smap = adapt_array_2d(squeeze(ref));
    smap=double(smap/max(abs(smap(:))));  

    nt=floor(ntviews/spokes_per_frame);
    clear kdatau dcfu traju
    for t = 1:nt
        kdatau(:,:,:,t) = kdata(:,(t-1)*spokes_per_frame+1:t*spokes_per_frame,:); % [nx,ntviews,ncoil] -> [nx,nspokes,ncoil,nt]
        dcfu(:,:,t) = dcf(:,(t-1)*spokes_per_frame+1:t*spokes_per_frame);
        traju(:,:,t) = traj(:,(t-1)*spokes_per_frame+1:t*spokes_per_frame);
    end
    
%     smap=double(smap/max(abs(smap(:))));
    param.E = MCNUFFT(traju,dcfu,smap);
    % Check w/ sqrt(dcf)
%     param.y=kdatau;
    param.y = double(kdatau.*repmat(reshape(sqrt(dcfu),[nx,spokes_per_frame,1,nt]),...
        [1,1,size(kdatau,3),1]));

    recon_nufft=param.E'*param.y;
    % parameters for reconstruction
    param.W = TV_Temp();
    param.lambda=option.lamb*max(abs(recon_nufft(:)));
    %param.lambda=0.0*max(abs(recon_nufft(:))); % iSENSE
    param.nite = option.niter;
    param.display=1;
    fprintf('\n GRASP reconstruction \n')
    tic
    recon_cs=recon_nufft;
    
    n = 1;
    bStop=0;
    error = 1;
    while (~bStop) 
%     for n=1:3
        
	    [recon_cs, recon_fval] = CSL1NlCg(recon_cs,param);
        cost_f_val(n) = recon_fval;
        if n>=2
            error = (cost_f_val(n)-cost_f_val(n-1))/cost_f_val(n-1);
        end
        if (abs(error)<2.5/100) & (n>3)
            bStop = 1;
        end
        n = n+1;
    end
    toc
    nufft_recon = recon_nufft;
    grasp_recon = recon_cs;
    b1 = [];
%     nufft_recon=flipdim(recon_nufft,1);
%     grasp_recon=flipdim(recon_cs,1);
else if strcmpi(option.method,'grog')
        nx = size(kdata,1);
        ntviews = size(kdata,2);
        Traj = traj;
%         Traj = Trajectory_GoldenAngle_GROG(ntviews,nx);
%         dcf = repmat(abs(linspace(-1,1,nx))',[1,ntviews]);
        nline = spokes_per_frame;

        if length(size(kdata))<4
            kdata = reshape(kdata,[size(kdata,1),size(kdata,2),1,size(kdata,3)]);
        end

        [Gx,Gy] = GROG.get_Gx_Gy(kdata,Traj);

        %Coil sensitivities estimation
        G = GROG.init(kdata,Traj,Gx,Gy,0);
        kref = GROG.interp(kdata,G,1);
        ref=squeeze(ifft2c_mri(kref));
        b1=adapt_array_2d(ref);
%         clear ref
        b1=single(b1/max(abs(b1(:))));
        %b1 = gpuArray(b1);
        
        %Calculate the DCF for nyquist sampling
        %Nqu=floor(bas/2*pi);
        G = GROG.init(kdata,Traj,Gx,Gy,0);
        DCF=reshape(G.weight,[sqrt(size(G.weight,1)),sqrt(size(G.weight,1))]);
        [nx,ntviews,nt,nc] = size(kdata);
        DCF=CropImg(DCF,nx,nx);
        
        %sort the data
        ntviews = size(kdata,2);
        nt=floor(ntviews/nline); %number of dynamic frames
        Traj=reshape(Traj(:,1:nt*nline),[nx,nline,nt]);
        kdata=reshape(kdata(:,1:nt*nline,:,:),[nx,nline,nt,nc]);
        [nx,ntviews,nt,nc] = size(kdata);
        
        %calculat weighting for iteration
        G = GROG.init(kdata,Traj,Gx,Gy,0);
        DCF_U=reshape(G.weight,[sqrt(size(G.weight,1)),sqrt(size(G.weight,1)),nt]);
        DCF_U=CropImg(DCF_U,nx,nx);
        DCF=repmat(DCF,[1,1,nt,nc]);
        DCF_U=repmat(DCF_U,[1,1,1,nc]);
        %Weighting=gpuArray(single(DCF./DCF_U));
        Weighting=DCF./DCF_U;
        
        %grog
        kdata = GROG.interp(kdata,G,1);
        mask=single(kdata~=0);
        kspace_cart = kdata;
        if option.GPU==1
            Weighting = gpuArray(Weighting);
            mask=gpuArray(mask);
            b1=gpuArray(b1);
            kdata = gpuArray(kdata);
        end
        
        bas=size(kdata,1)/2;
        
        param.E=Emat_GROG2D(mask,b1,single(Weighting));
        %keyboard;
        
        %clear mask b1 DCF DCF_U kref 
        
        param.y=single(kdata).*sqrt(single(Weighting));
        
        recon_cs=param.E'*param.y;
        data=(single(gather(recon_cs/max(recon_cs(:)))));

        param.TV=TV_Temp;
        Weight1=option.lamb;

        param.TVWeight=max(abs(recon_cs(:)))*Weight1;
        param.nite = option.niter;
        param.display = 1;
        n = 1;
        bStop=0;
        error = 1;
        while ~bStop
    %         for n=1:3
            [recon_cs,recon_fval] = CSL1NlCg_GROG(recon_cs,param);
            cost_f_val(n) = recon_fval;
            if n>=2
                error = (cost_f_val(n)-cost_f_val(n-1))/cost_f_val(n-1);
            end
            if abs(error)<2.5/100
                bStop = 1;
            end
            n = n+1;
        end
%         data(:,:,:,2)=gather(single(recon_cs/max(recon_cs(:))));
        data(:,:,:,2)=gather(single(recon_cs));
        data=CropImg(data,bas,bas);
        nufft_recon = data(:,:,:,1);
        grasp_recon = data(:,:,:,2);
        b1 = gather(CropImg(b1,bas,bas));
end
end


