% chin feb 5 2013

% emd fault identification with real data

clear
clc
addpath Utils
addpath emd_flow

% set parameters k, B, thrsh

load seismic_fault_real

k = 15;
B = 150;
shft = 5;
offset1 = 125; offset2 = 20;
Zfaultn = a(21:95,126:200); Zfaultn = Zfaultn/max(abs(Zfaultn(:)));
Zfaultn = double(Zfaultn);

% k = 15;
% B = 150;
% shft = 5;
% offset1 = 125; offset2 = 125;
% Zfaultn = a((offset2+1):250,(offset1+1):200); Zfaultn = Zfaultn/max(abs(Zfaultn(:)));
% Zfaultn = double(Zfaultn);


figure(10), clf
subplot(1,2,1), imagesc(Zfaultn), axis image %, caxis([-1 1])
subplot(1,2,2), imagesc(Zfaultn), axis image %, caxis([-1 1])
hold on
scatter(faults2(:,1)-offset1,faults2(:,3)-offset2,50,[0 0 0],'filled')

%return
opts.verbose = true;
%
mags = double(Zfaultn.^2);
supp = emd_flow(mags,k,B,opts);
supp = double(supp);

%return

thrshvec = 1.5;

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
idxoutlier = find(Dsort > 2*median(Dsort));
idxkeep = setdiff(1:length(idx),idxoutlier);
v1sort_keep = v1sort(idxkeep);
v2sort_keep = v2sort(idxkeep);

figure(2), clf, 
subplot(1,3,1),
imagesc(Zfaultn), axis image
axisfortex('','Noisy input','')
rmaxis
subplot(1,3,2),
imagesc(Zfaultn), axis image
hold on
scatter(faults2(:,1)-offset1,faults2(:,3)-offset2,50,[0 0 0],'filled')
axisfortex('','Human labels','')
rmaxis
subplot(1,3,3),
imagesc(Zfaultn), axis image
hold on
scatter(v2sort_keep,v1sort_keep,50,'kd','filled')
axisfortex('','Automatic','')
rmaxis

disp(thrsh)
pause(2)

end