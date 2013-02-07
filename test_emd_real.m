% chin feb 5 2013

% emd fault identification with real data

clear
clc
addpath Utils
addpath emd_flow

% set parameters k, B, thrsh

% k = 15, B = 400 for a test case that breaks monotonicity?
k = 15;
B = 400;
shft = 5;
thrsh = 0.5*shft;

load seismic_fault_real

Zfaultn = a(:,51:end); Zfaultn = Zfaultn/max(abs(Zfaultn(:)));

Zfaultn = sign(Zfaultn).*(Zfaultn.^2);

figure(10), clf
subplot(1,2,1), imagesc(Zfaultn), axis image %, caxis([-1 1])
subplot(1,2,2), imagesc(Zfaultn), axis image %, caxis([-1 1])
hold on
scatter(faults2(:,1)-50,faults2(:,3),20,'filled')

%return


mags = double(Zfaultn.^2);
supp = emd_flow(mags,k,B,true);
supp = double(supp);

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
subplot(1,2,1),
imagesc(Zfaultn), axis image
axisfortex('','Noisy input','')
rmaxis
subplot(1,3,3),
imagesc(Zfaultn), axis image
hold on
plot(vy,vx,'k--','LineWidth',4)
axisfortex('','Labeled Fault(s)','')
rmaxis
