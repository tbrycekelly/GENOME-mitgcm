% ------------------------- %
% Read MITgcm input files
% Written By: Taylor Shropshire 
% Original Date: 7/5/18
% --------------------------- %
clear;clc;close all

% Params
path1='/nexsan/people/taylor/projects/mit_gcm/MITgcm2/MITgcm/verification/gom_test_sim/input/flow_fields/';
xdim=104;
ydim=72;
zdim=29;

% Get file names 
fnames_uvel=dir([path1,'uvel/*.data']);
fnames_vvel=dir([path1,'vvel/*.data']);
fnames_wvel=dir([path1,'wvel/*.data']);
fnames_tmp=dir([path1,'tmp/*.data']);
fnames_sal=dir([path1,'sal/*.data']);

% Read each daily average flow field
uvel_all=[];vvel_all=[];wvel_all=[];
tmp_all=[];sal_all=[];

%for tt=1:length(fnames_uvel)
tt=1;
display(['Day = ', num2str(tt)]);
tic
file1=[path1,'uvel/',fnames_uvel(tt).name];
file2=[path1,'vvel/',fnames_vvel(tt).name];
file3=[path1,'wvel/',fnames_wvel(tt).name];
file4=[path1,'tmp/',fnames_tmp(tt).name];
file5=[path1,'sal/',fnames_sal(tt).name];

fid1=fopen(file1,'r','b');uvel=fread(fid1,'float');fclose(fid1);
fid2=fopen(file2,'r','b');vvel=fread(fid2,'float');fclose(fid2);
fid3=fopen(file3,'r','b');wvel=fread(fid3,'float');fclose(fid3);
fid4=fopen(file4,'r','b');tmp=fread(fid4,'float');fclose(fid4);
fid5=fopen(file5,'r','b');sal=fread(fid5,'float');fclose(fid5);

% Read each line of binary file
uvel2=[];vvel2=[];wvel2=[];tmp2=[];sal2=[];
ind1=1;
ind2=xdim;
for zz=1:zdim
for xx=1:ydim
uvel2(:,xx,zz)=uvel(ind1:ind2);
vvel2(:,xx,zz)=vvel(ind1:ind2);
wvel2(:,xx,zz)=wvel(ind1:ind2);
tmp2(:,xx,zz)=tmp(ind1:ind2);
sal2(:,xx,zz)=sal(ind1:ind2);
ind1=ind1+xdim;
ind2=ind2+xdim;
end
end
uvel=uvel2;clear uvel2;uvel(find(uvel==0))=NaN;
vvel=vvel2; clear vvel2;vvel(find(vvel==0))=NaN;
wvel=wvel2;clear wvle2;wvel(find(wvel==0))=NaN;
tmp=tmp2; clear tmp2;tmp(find(tmp==0))=NaN;
sal=sal2; clear sal2;sal(find(sal==0))=NaN;

% Rotate 
uvel2=[];vvel2=[];wvel2=[];tmp2=[];sal2=[];
for zz=1:size(uvel,3)
uvel2(:,:,zz)=uvel(:,:,zz)';
vvel2(:,:,zz)=vvel(:,:,zz)';
wvel2(:,:,zz)=wvel(:,:,zz)';
tmp2(:,:,zz)=tmp(:,:,zz)';
sal2(:,:,zz)=sal(:,:,zz)';
end 
uvel=uvel2;clear uvel2;vvel=vvel2; clear vvel2;
wvel=wvel2;clear wvle2;tmp=tmp2; clear tmp2;sal=sal2; clear sal2

uvel_all(:,:,:,tt)=uvel;
vvel_all(:,:,:,tt)=uvel;
wvel_all(:,:,:,tt)=uvel;
tmp_all(:,:,:,tt)=uvel;
sal_all(:,:,:,tt)=uvel;
toc 
%end 

uvel_field=mean(uvel_all,4);
vvel_field=mean(vvel_all,4);
wvel_field=mean(wvel_all,4);
tmp_field=mean(tmp_all,4);
sal_field=mean(sal_all,4);

figure(1)
pcolor(uvel_field(:,:,1));shading flat;colorbar;caxis([-1 1])
title('u-vel (m^-^3 /s)')



