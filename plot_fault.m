% function to visualize labeled 3D fault profiles

% chin may 24 2013

% input: point cloud (in the form of Nx3 array), parameters

% assume for now that the data is a cube

function probdensity = plot_fault(point_cloud, opt)


[X,Y,Z] = meshgrid(1:opt.datalim);
probdirac = 0*X;
for ii=1:size(point_cloud,1)
  probdirac(point_cloud(ii,1),point_cloud(ii,2),point_cloud(ii,3)) = 1;
end

% smooth the prob density

opt.amp = 1;
opt.sigmasq = 0.5*opt.kernelsize;
gaussKernel = zeros(opt.kernelsize,opt.kernelsize,opt.kernelsize);
mid = ceil(0.5*opt.kernelsize);
for x = 1:opt.kernelsize, 
    for y=1:opt.kernelsize, 
        for z=1:opt.kernelsize
radiusSquared = (x-mid)^2 + (y-mid)^2 + (z-mid)^2;
gaussKernel(x, y, z) = opt.amp * exp(-radiusSquared/opt.sigmasq);
        end; 
    end; 
end
probdensity = imfilter(probdirac,gaussKernel);

figure; hold on;
for a = 10.^(1:4) % 'a' defines the isosurface limits
p = patch(isosurface(X,Y,Z,probdensity,max(max(max(probdensity)))/a)); % isosurfaces at max(V)/a
isonormals(X,Y,Z,probdensity,p); % plot the surfaces
set(p,'FaceColor','red','EdgeColor','none'); % set colors
end
alpha(.1); % set the transparency for the isosurfaces
xlim([0 opt.datalim]); ylim([0 opt.datalim]); zlim([0 opt.datalim]);
axis square
box on
%axisfortex('','','')
