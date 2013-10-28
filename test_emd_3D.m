% chin may 23 2013
% emd fault identification with real data

%clear
clc
addpath Utils
addpath emd_flow-master/

% set parameters k, B, thrsh

load faults_image
k = 15;
B = 150;
offset1 = 125; offset2 = 20;
idx = 31:70;

human_faults = [];
auto_faults = [];

opts.outdegree_vertical_distance = 10;

for kk=idx

inds2 = find(faults_ind_small(:,2)==kk);
faults2 = faults_ind_small(inds2,:);
a = squeeze(seis_image(1:250,1:200,kk));
Zfaultn = a(21:95,126:200); 
Zfaultn = Zfaultn/max(abs(Zfaultn(:)));
Zfaultn = double(Zfaultn);

%return
opts.verbose = false;
mags = double(Zfaultn.^2);
supp = emd_flow(mags,k,B,opts);
supp = double(supp);

%return
thrshvec = 1.5;
medfact = 1.75;

for thrsh = thrshvec

% label the faults
[i,j] = find(supp); % determine the flows
i = reshape(i,k,[]);
j = reshape(j,k,[]);
flowdiff = abs(diff(i,1,2));
m = median(flowdiff,2);
[i1,i2] = find(flowdiff > repmat(m,1,size(flowdiff,2)) + thrsh);

v1 = diag(i(i1,i2));
v2 = diag(j(i1,i2));
[v1sort, idx] = sort(v1,'descend');
v2sort = v2(idx);

% some basic outlier rejection
% if there's a solitary marked point somewhere, 
% discard it

X = [v1sort v2sort]';
D = L2_distance(X,X);
D = D + max(D(:))*diag(ones(length(D),1));
[Dsort,idx] = min(D,[],2);
idxoutlier = find(Dsort > medfact*median(Dsort));
idxkeep = setdiff(1:length(idx),idxoutlier);
v1sort_keep = v1sort(idxkeep);
v2sort_keep = v2sort(idxkeep);

figure(2), clf, 
subplot(1,3,1),
imagesc(Zfaultn), axis image, caxis([-0.7 0.7])
axisfortex('','Noisy input','')
rmaxis
subplot(1,3,2),
imagesc(Zfaultn), axis image, caxis([-0.7 0.7])
hold on
scatter(faults2(:,1)-offset1,faults2(:,3)-offset2,50,[0 0 0],'filled')
axisfortex('','Human labels','')
rmaxis
subplot(1,3,3),
imagesc(Zfaultn), axis image, caxis([-0.7 0.7])
hold on
scatter(v2sort_keep,v1sort_keep,50,'kd','filled')
axisfortex('','Automatic','')
rmaxis

% build 3D human-label fault profile
a = [faults2(:,1)-offset1,faults2(:,3)-offset2];
c = a(find((a(:,1)>0).*(a(:,1)<76).*(a(:,2)>0).*(a(:,2)<76)),:);
human_faults = [human_faults; c,repmat(kk,size(c,1),1)];

% build 3D automatic fault profile
c = [v2sort_keep,v1sort_keep];
auto_faults = [auto_faults; c,repmat(kk,size(c,1),1)];

disp(kk)
pause(0.05)

end



end

opt.filt = 1;
opt.filtfact = 1;
%%%%% 3D outlier filtering
%%%%% exploit cons

if opt.filt == 1
X = auto_faults';
D = L2_distance(X,X);
D = D + max(D(:))*diag(ones(length(D),1));

[Dsort,idx] = min(D,[],2);
idxoutlier = find(Dsort > medfact*median(Dsort));
idxkeep = setdiff(1:length(idx),idxoutlier);

Dsort = sort(D);

idxoutlier1 = find(Dsort(1,:) > opt.filtfact*median(Dsort(1,:)));
idxoutlier2 = find(Dsort(2,:) > opt.filtfact*median(Dsort(2,:)));
idxoutlier3 = find(Dsort(3,:) > opt.filtfact*median(Dsort(3,:)));
%idxoutlier3 = 1:size(Dsort,2);

idxoutlier = intersect(intersect(idxoutlier2,idxoutlier3),idxoutlier1);

idxkeep = setdiff(1:size(D,2),idxoutlier);

auto_faults = auto_faults(idxkeep,:);
end







figure(30), clf
%subplot(1,2,1)
scatter3(human_faults(:,3),(human_faults(:,1)),human_faults(:,2),25,'fill')
%axis([0 76 0 76 0 76])
view([-45 30])
xlim([26 75]); ylim([0 76]); zlim([0 76]);
axis square
box on
axisfortex('Human Labels','X','Y')
%subplot(1,2,2)
figure(31), clf
scatter3(auto_faults(:,3),(auto_faults(:,1)),auto_faults(:,2),25,'fill')
view([-45 30])
xlim([26 75]); ylim([0 76]); zlim([0 76]);
axis square
box on
axisfortex('Auto Labels','X','Y')


%%%%%%%%%%%%%%%%%%%%%%%%
% plane fitting
%%%%%%%%%%%%%%%%%%%%%%%%


human_faults_aug = [human_faults ones(size(human_faults,1),1)];
[u,s,v] = svd(human_faults_aug,'econ');
P = v(:,4);
pos = random('unid',75,[10000 2]);
Z = (-1/P(3))*(pos*P(1:2) + P(4));
pos3 = [pos Z];
figure(32), clf 
scatter3(pos3(:,3),(pos3(:,1)),pos3(:,2),25,'fill')
view([-45 30])
xlim([26 75]); ylim([0 76]); zlim([0 76]);
axis square
box on
axisfortex('Plane Fit (Human)','X','Y')


auto_faults_aug = [auto_faults ones(size(auto_faults,1),1)];
[u,s,v] = svd(auto_faults_aug,'econ');
P = v(:,4);
pos = random('unid',75,[10000 2]);
Z = (-1/P(3))*(pos*P(1:2) + P(4));
pos3 = [pos Z];
figure(33), clf 
scatter3(pos3(:,3),(pos3(:,1)),pos3(:,2),25,'fill')
view([-45 30])
xlim([26 75]); ylim([0 76]); zlim([0 76]);
axis square
box on
axisfortex('Plane Fit (Auto)','X','Y')
