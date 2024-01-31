function [simImg,mask,parMap,temp,ID,aif,T10,aifci_1s,cts] = gen_simIMG(data,S0)


% idx = round(length(data)*rand());
idx = 3;

par_var = 0.1;
par_var_t = 0.2;

t = linspace(0,150,22);
%t = linspace(0,15,132);
mask = data{idx}.mask;

aif = data{idx}.AIF;
% aif_idx = round(length(data)*rand());
% aif = data{aif_idx}.AIF;

if ~exist('S0')
% S0 = abs(data{idx}.S0);
S0 = data{idx}.S0;
end
ID = data{idx}.ID;
smap = double(data{idx}.smap);
disp(ID);

if ~isfield(mask,'heart')
    mask.heart = logical(zeros(size(S0,1),size(S0,2)));
end
if ~isfield(mask,'liver')
    mask.liver = logical(zeros(size(S0,1),size(S0,2)));
end
if ~isfield(mask,'skin')
    mask.skin = logical(zeros(size(S0,1),size(S0,2)));
end
if ~isfield(mask,'muscle')
    mask.muscle = logical(zeros(size(S0,1),size(S0,2)));
end
if ~isfield(mask,'benign')
    mask.benign = logical(zeros(size(S0,1),size(S0,2)));
end
if ~isfield(mask,'malignant')
    mask.malignant = logical(zeros(size(S0,1),size(S0,2)));
end
if ~isfield(mask,'heart_blood')
    mask.heart_blood = mask.heart;
end


nbase = find((aif<0.15)&(t<300));
nbase = nbase(end);


t_1s = 0:t(end);
aifci = interp1(t,aif,t_1s,'pchip');

[nx,ny] = size(S0);
%logIdx = gen_downSample_logIdx(t,ti);


T1.glandular = 1.324;
T1.malignant = 1.5;
T1.benign = 1.4;
T1.liver = 0.81;
T1.heart = 1.68;
T1.muscle = 1.41;
T1.skin = 0.85; %Check
% T1.skin = 1;
T1.vasc = 1.93;

T10 = zeros(nx,ny);
T10(mask.glandular) = T1.glandular;
T10(mask.malignant) = T1.malignant;
T10(mask.benign) = T1.benign;
T10(mask.liver) = T1.liver;
T10(mask.heart) = T1.heart;
T10(mask.muscle) = T1.muscle;
T10(mask.skin) = T1.skin;
T10(T10==0) = 1e-8;

% [Ve,Vp,Fp,PS]
% Ve: 20~30% (.2~.3)
% Vp: 2~10% (.02~.1)
% Fp: 0.5~2.5 (min-1)
% PS: 0.01~0.8 (min-1)


% p0.glandular = [0.0280,0.161;0.013,0.041;0.034/60,0.132/60;0.236/60,1.269/60];
% p0.malignant = [0.024,0.270;0.013,0.241;0.290/60,2.47/60;0.206/60,1.287/60];
% p0.benign = [0.052,0.258;0.017,0.132;0.125/60,0.319/60;0.172/60,1.272/60];
% p0.muscle = [0.058,0.169;0.011,0.05;0.070/60,0.284/60;0.141/60,1.253/60];
% p0.skin = [0.06,0.145;0.02,0.025;0.365/60,1.563/60;0.069/60,0.148/60]; % Skin vp<0.01 => nan
% p0.heart = [0.190,0.369;0.127,0.239;3/60,3/60;0.447/60,0.870/60];
% p0.liver = [0.01,0.014;0.094,0.240;0.179/60,0.352/60;0.221/60,1.282/60];

p0.glandular = [0.049,0.162;0.01,0.026;0.061/60,0.121/60;0.092/60,1.272/60];
p0.malignant = [0.010,0.263;0.103,0.317;0.274/60,0.905/60;0.130/60,.464/60];
p0.benign = [0.052,0.262;0.017,0.135;0.142/60,0.356/60;0.191/60,1.289/60];
p0.muscle = [0.080,0.155;0.011,0.02;0.075/60,0.243/60;0.118/60,1.222/60];
% p0.skin = [0.062,0.151;0.005,0.01;0.285/60,0.743/60;0.074/60,0.154/60]; % Skin vp<0.01 => nan
p0.skin = [0.062,0.151;0.01,0.012;0.285/60,0.743/60;0.074/60,0.154/60]; % Skin vp<0.01 => nan
p0.heart = [0.189,0.329;0.146,0.317;2/60,3/60;0.439/60,0.909/60];
p0.liver = [0.01,0.017;0.122,0.386;0.170/60,0.356/60;0.131/60,0.989/60];


gland_ktrans = [0.034,0.094]./60;
malig_ktrans = [0.110,0.231]./60;
benign_ktrans = [0.102,0.214]./60;
muscle_ktrans = [0.066,0.106]./60;
skin_ktrans = [0.07,0.120]./60;
liver_ktrans = [0.110,0.199]./60;
heart_ktrans = [0.408,0.769]./60;

fun_ktrans = @(x) x(3).*(1-exp(-x(4)./x(3)));
p0_glandular = p0.glandular(:,1)+(p0.glandular(:,2)-p0.glandular(:,1)).*rand(4,1);
par_ktrans = fun_ktrans(p0_glandular);
while ((par_ktrans)<gland_ktrans(1)||(par_ktrans)>gland_ktrans(2))
    p0_glandular = p0.glandular(:,1)+(p0.glandular(:,2)-p0.glandular(:,1)).*rand(4,1);
    par_ktrans = fun_ktrans(p0_glandular);
end
p0_malignant = p0.malignant(:,1)+(p0.malignant(:,2)-p0.malignant(:,1)).*rand(4,1);
par_ktrans = fun_ktrans(p0_malignant);
while ((par_ktrans)<malig_ktrans(1)||(par_ktrans)>malig_ktrans(2))
    p0_malignant = p0.malignant(:,1)+(p0.malignant(:,2)-p0.malignant(:,1)).*rand(4,1);
    par_ktrans = fun_ktrans(p0_malignant);
end
p0_benign = p0.benign(:,1)+(p0.benign(:,2)-p0.benign(:,1)).*rand(4,1);
par_ktrans = fun_ktrans(p0_benign);
while ((par_ktrans)<benign_ktrans(1)||(par_ktrans)>benign_ktrans(2))
    p0_benign = p0.benign(:,1)+(p0.benign(:,2)-p0.benign(:,1)).*rand(4,1);
    par_ktrans = fun_ktrans(p0_benign);
end
p0_muscle = p0.muscle(:,1)+(p0.muscle(:,2)-p0.muscle(:,1)).*rand(4,1);
par_ktrans = fun_ktrans(p0_muscle);
while ((par_ktrans)<muscle_ktrans(1)||(par_ktrans)>muscle_ktrans(2))
    p0_muscle = p0.muscle(:,1)+(p0.muscle(:,2)-p0.muscle(:,1)).*rand(4,1);
    par_ktrans = fun_ktrans(p0_muscle);
end
p0_skin = p0.skin(:,1)+(p0.skin(:,2)-p0.skin(:,1)).*rand(4,1);
par_ktrans = fun_ktrans(p0_skin);
while ((par_ktrans)<skin_ktrans(1)||(par_ktrans)>skin_ktrans(2))
    p0_skin = p0.skin(:,1)+(p0.skin(:,2)-p0.skin(:,1)).*rand(4,1);
    par_ktrans = fun_ktrans(p0_skin);
end
p0_heart = p0.heart(:,1)+(p0.heart(:,2)-p0.heart(:,1)).*rand(4,1);
par_ktrans = fun_ktrans(p0_heart);
while ((par_ktrans)<heart_ktrans(1)||(par_ktrans)>heart_ktrans(2))
    p0_heart = p0.heart(:,1)+(p0.heart(:,2)-p0.heart(:,1)).*rand(4,1);
    par_ktrans = fun_ktrans(p0_heart);
end
p0_liver = p0.liver(:,1)+(p0.liver(:,2)-p0.liver(:,1)).*rand(4,1);
par_ktrans = fun_ktrans(p0_liver);
while ((par_ktrans)<liver_ktrans(1)||(par_ktrans)>liver_ktrans(2))
    p0_liver = p0.liver(:,1)+(p0.liver(:,2)-p0.liver(:,1)).*rand(4,1);
    par_ktrans = fun_ktrans(p0_liver);
end
parMap = zeros(nx,ny,4);

for i = 1:4
    temp = zeros(nx,ny);
    randMap = p0_glandular(i)*((1-par_var)+(par_var*2)*rand(nx,ny));
    temp(mask.glandular) = randMap(mask.glandular);
    
    randMap = p0_malignant(i)*((1-par_var_t)+(par_var_t*2)*rand(nx,ny));
    temp(mask.malignant) = randMap(mask.malignant);
    
    randMap = p0_benign(i)*((1-par_var_t)+(par_var_t*2)*rand(nx,ny));
    temp(mask.benign) = randMap(mask.benign);
    
    %temp(mask.fat) = p0.fat(i);
    randMap = p0_liver(i)*((1-par_var)+(par_var*2)*rand(nx,ny));
    temp(mask.liver) = randMap(mask.liver);
    
    %temp(mask.myocard) = p0.myocard(i);
    randMap = p0_muscle(i)*((1-par_var)+(par_var*2)*rand(nx,ny));
    temp(mask.muscle) = randMap(mask.muscle);
    
    randMap = p0_skin(i)*((1-par_var)+(par_var*2)*rand(nx,ny));
    temp(mask.skin) = randMap(mask.skin);
    %temp = temp.*(0.8+(0.4*rand(nx,ny)));

    randMap = p0_heart(i)*((1-par_var)+(par_var*2)*rand(nx,ny));
    temp(mask.heart_blood) = randMap(mask.heart_blood);
    parMap(:,:,i) = temp;
end


aifci_1s = zeros(nx,ny,length(t_1s));
[rIdx,cIdx] = find(mask.heart_blood==1); %4s delay
for i = 1:length(rIdx)
    aifci_1s(rIdx(i),cIdx(i),:) = aifci([1:4,1:end-4]);
end
[rIdx,cIdx] = find(mask.liver==1); %8s delay
for i = 1:length(rIdx)
    aifci_1s(rIdx(i),cIdx(i),:) = aifci([1:8,1:end-8]);
end

[rIdx,cIdx] = find(mask.glandular==1);%15s delay
for i = 1:length(rIdx)
    aifci_1s(rIdx(i),cIdx(i),:) = aifci([1:15,1:end-15]);
end

[rIdx,cIdx] = find(mask.malignant==1); %8s delay
for i = 1:length(rIdx)
    aifci_1s(rIdx(i),cIdx(i),:) = aifci([1:8,1:end-8]);
end
[rIdx,cIdx] = find(mask.benign==1); %13s delay
for i = 1:length(rIdx)
    aifci_1s(rIdx(i),cIdx(i),:) = aifci([1:13,1:end-13]);
end

[rIdx,cIdx] = find(mask.muscle==1); %10s delay
for i = 1:length(rIdx)
    aifci_1s(rIdx(i),cIdx(i),:) = aifci([1:10,1:end-10]);
end
[rIdx,cIdx] = find(mask.skin==1); %10s delay
for i = 1:length(rIdx)
    aifci_1s(rIdx(i),cIdx(i),:) = aifci([1:10,1:end-10]);
end
% bgd_mask = parMap(:,:,4)==0;


ti = 0:0.1:t(end);

logIdx = zeros(1,length(ti));
start_idx = 1;
for i = 1:length(t)
    for j = start_idx:length(ti)
        if t(i)<=ti(j)
            logIdx(j) = 1;
            start_idx = j;
            break
        end
    end
end
logIdx = logical(logIdx);

aifci_Map = zeros(nx,ny,length(ti));
[rIdx,cIdx] = find(parMap(:,:,1)>0);
for i = 1:length(rIdx)
    temp_aif = squeeze(aifci_1s(rIdx(i),cIdx(i),:));
    aifci_Map(rIdx(i),cIdx(i),:) = interp1(t_1s,temp_aif,ti,'pchip');
end


parMap(parMap==0) = 1e-8;

for i = 1:4
%     parMap(:,:,i) = imgaussfilt(parMap(:,:,i),1.7);
    parMap(:,:,i) = imgaussfilt(parMap(:,:,i),0.1);
end


ve = parMap(:,:,1);
vp = parMap(:,:,2);
fp = parMap(:,:,3);
ktrans = parMap(:,:,4);


Ce = zeros(size(parMap,1),size(parMap,2),size(aifci_Map,3));
Cp = zeros(size(parMap,1),size(parMap,2),size(aifci_Map,3));

for i = 2:size(aifci_Map,3)
    dt = ti(i)-ti(i-1);
    dcp = fp.*aifci_Map(:,:,i-1) - (fp+ktrans).*Cp(:,:,i-1) + ktrans.*Ce(:,:,i-1);
    dce = ktrans.*Cp(:,:,i-1)-ktrans.*Ce(:,:,i-1);
    Cp(:,:,i) = Cp(:,:,i-1) + dcp.*dt./vp;
    Ce(:,:,i) = Ce(:,:,i-1) + dce.*dt./ve;
end

cts = Cp.*repmat(vp,[1,1,length(ti)]) + Ce.*repmat(ve,[1,1,length(ti)]);
cts = cts(:,:,logIdx);

nan_mask = sum(isnan(cts),3);
[rIdx,cIdx] = find(nan_mask);
disp(['Found ',num2str(length(rIdx)),' NaN voxels']);
for i = 1:length(rIdx)
    cts(rIdx(i),cIdx(i),:) = 0;
end

TR = 4.87e-3;
theta = 10 * pi /180;
r1 = 4.3;
T1_t = 1./(1./repmat(T10,[1,1,size(cts,3)]) + r1*cts);

Eh = (1 - exp(-TR./T1_t)).*sin(theta)./ (1-(exp(-TR./T1_t).*cos(theta)));
Eh = Eh ./ repmat(mean(Eh(:,:,1:nbase),3),[1,1,size(Eh,3)]);
%Eh(mask_enh) = 1e-8;


simImg = Eh .* repmat(S0,[1,1,size(Eh,3)]);
% aif_map = aifci_Map(:,:,logIdx);

nan_mask = sum(isnan(simImg),3);
[rIdx,cIdx] = find(nan_mask);
disp(['Sig-Found ',num2str(length(rIdx)),' NaN voxels']);
if nnz(nan_mask)>0
    close all;
    figure;
    imshow(abs(S0),[]);colormap(gray);
    hold on;freezeColors();
    h = imshow(nan_mask>0,[]);colormap(jet);
    set(h,'AlphaData',nan_mask>0);

    keyboard;
end

%simImg = simImg(:,:,logIdx);
%fsimImg = simImg / max(simImg(:));