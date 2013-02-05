% chin jan 29 2013

% automatic fault detection

% this file makes test images with a single (slanted) fault
function [Z,Zfault,Zfaultn] = make_test_image(n,w,sg,shft,theta,flg)

%n = 150; shft = 5;
[X,Y] = meshgrid(-w:w,-n:n);
if flg==0
    Z = zeros(2*n+1,2*w+1);
    Z(round(0.3*n):round(0.3*n),:) = 1;
    Z(round(0.6*n):round(0.6*n),:) = 1;
    Z(round(0.8*n):round(0.8*n),:) = -1;
    Z(round(1.5*n):round(1.5*n),:) = 1;
    Z(round(1.7*n):round(1.7*n),:) = -1;
elseif flg ==1
      Znew = load('wavy_lines.mat');  
      Z = Znew.x;
      Z(50:end,:) = -1*Z(50:end,:);
elseif flg ==2
      Znew = load('wavy_lines2.mat');  
      Z = Znew.x;
      Z(50:end,:) = -1*Z(50:end,:);
end

%theta = -75; 
m = tan(theta*pi/180);
Zbit = double((Y - m*X)>0);
%[I,J] = find(Zbit > 0);
Zshft = 0*Z; Zshft(1:(2*n-shft),:) = Z((shft+1):2*n,:);
Zfault = (1-Zbit).*Z + Zbit.*Zshft;
Zfaultn = Zfault + sg*randn(size(Zfault));

% figure(1), clf
% subplot(1,3,1), imagesc(Z), colormap(gray)
% subplot(1,3,2), imagesc(Zfaultn), colormap(gray)
% subplot(1,3,3), imagesc(Fx)
