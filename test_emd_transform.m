% chin jan 29 2013

clear
clc
addpath Utils
addpath /Users/Chin/Documents/Chin/Acads/modelcs/modelcs_debug/EMD/
addpath /Users/Chin/Documents/Chin/Acads/modelcs/modelcs_debug/EMD/emd_flow/

sg = 0.7;
n = 50;
w = 25;
shft = 3;
theta = 77;

flg = 0;

[Z,Zfault,Zfaultn] = make_test_image(n,w,sg,shft,theta,flg);
% calculate EMD model projection
if flg==0
k = 5; 
thrsh = 0.5*shft;
B = round(1.1*k*shft);
elseif flg == 1
   k = 3;
   B = 160;
   thrsh = 0.5*shft;
elseif flg == 2
    k = 4;
    B = 200+shft;
    thrsh = 0.5*shft;
end

mags = Zfaultn.^2;
supp = emd_flow(mags,k,B,true);
supp = double(supp);
% 
% figure(1), clf 
% subplot(1,3,1)
% imagesc(Zfaultn)
% subplot(1,3,2)
% imagesc(supp)

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
vx = [v1sort(1:end-1)'; v1sort(2:end)'];
vy = [v2sort(1:end-1)'; v2sort(2:end)'];

figure(2), clf, 
subplot(1,3,1),
imagesc(Zfault), axis image
axisfortex('','Original','')
rmaxis
subplot(1,3,2),
imagesc(Zfaultn), axis image
axisfortex('','Noisy input','')
rmaxis
subplot(1,3,3),
imagesc(Zfaultn), axis image
hold on
plot(vy,vx,'k--','LineWidth',4)
axisfortex('','Labeled Fault','')
rmaxis
