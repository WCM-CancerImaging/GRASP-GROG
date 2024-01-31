% Reshaping k-space

[nx,ntview,nc] = size(kdata);
kdata = reshape(kdata,[nx,ntview,1,nc]);

%%
traj = Trajectory_GoldenAngle_GROG(ntview,nx); %for grog
% traj = Trajectory_GoldenAngle_GROG(ntview,nx)./nx; %for nufft
dcf = repmat(abs(linspace(-nx/2,nx/2,nx))',[1,ntview]);

spokes_per_frame = 13;
option.method = 'grog';
option.GPU = 1;
option.lamb = 0.01;
option.niter = 5;
grasp = run_iter_grasp(kdata,traj,dcf,spokes_per_frame,option);

