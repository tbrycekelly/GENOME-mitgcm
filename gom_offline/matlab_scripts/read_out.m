% ------------------------- %
% Read NetCDF MITgcm Simulation Output 
% Written By: Taylor Shropshire 
% Original Date: 4/4/17
% --------------------------- %
clear;clc;close all

path_in1='/nexsan/people/taylor/projects/mit_gcm/MITgcm2/MITgcm/verification/gom_test_sim/expt001/';
dir_names=dir([path_in1,'*tile*']);

% Get time dim 
dd=1;pp=1;
ptr_files=dir([path_in1,dir_names(dd).name,'/*ptracers*']);
nc_file1=[path_in1,dir_names(dd).name,'/',ptr_files(pp).name];
ncid1=netcdf.open(nc_file1,'NOWRITE');
varid=netcdf.inqVarID(ncid1,'T');td=netcdf.getVar(ncid1,varid);td=double(td);td=length(td);
netcdf.close(ncid1);

t_start=1;              
num_days=td-1; 
ptr_file_start=1;

% **** Read NetCDF state variable files ****
NO_all=[];NH_all=[];SI_all=[];DON_all=[];PON_all=[];OP_all=[];
sp_all=[];lp_all=[];sz_all=[];lz_all=[];pz_all=[];
chl_all=[];tmp_all=[];l_lim_sp_all=[];l_lim_lp_all=[];
xc_all2=[];yc_all2=[];

for dd=1:length(dir_names)
tic
ptr_files=dir([path_in1,dir_names(dd).name,'/*ptracers*']);
diag_files=dir([path_in1,dir_names(dd).name,'/*bio_diagnostics*']);
grid_file=dir([path_in1,dir_names(dd).name,'/*grid*']);

for pp=ptr_file_start:length(ptr_files)
display(['Reading File : ',dir_names(dd).name,'/',ptr_files(pp).name])

% Open netCDF file
nc_file1=[path_in1,dir_names(dd).name,'/',ptr_files(pp).name];
nc_file2=[path_in1,dir_names(dd).name,'/',diag_files(pp).name];
nc_file3=[path_in1,dir_names(dd).name,'/',grid_file(1).name];
ncid1=netcdf.open(nc_file1,'NOWRITE');
ncid2=netcdf.open(nc_file2,'NOWRITE');
ncid3=netcdf.open(nc_file3,'NOWRITE');

% Get field dimensions
varid=netcdf.inqVarID(ncid1,'T');td=netcdf.getVar(ncid1,varid);td=double(td);
varid=netcdf.inqVarID(ncid3,'drF');zd=netcdf.getVar(ncid3,varid);zd=double(zd);
varid=netcdf.inqVarID(ncid3,'XC');xc=netcdf.getVar(ncid3,varid);xc=double(xc);
varid=netcdf.inqVarID(ncid3,'YC');yc=netcdf.getVar(ncid3,varid);yc=double(yc);
xd=size(xc,1);yd=size(xc,2);td=length(td);zd=length(zd);

% Get state variables
varid=netcdf.inqVarID(ncid1,'NO');NO=netcdf.getVar(ncid1,varid,[0,0,0,t_start],[xd,yd,zd,num_days]);NO=double(NO);
varid=netcdf.inqVarID(ncid1,'NH');NH=netcdf.getVar(ncid1,varid,[0,0,0,t_start],[xd,yd,zd,num_days]);NH=double(NH);
varid=netcdf.inqVarID(ncid1,'SI');SI=netcdf.getVar(ncid1,varid,[0,0,0,t_start],[xd,yd,zd,num_days]);SI=double(SI);
varid=netcdf.inqVarID(ncid1,'DON');DON=netcdf.getVar(ncid1,varid,[0,0,0,t_start],[xd,yd,zd,num_days]);DON=double(DON);
varid=netcdf.inqVarID(ncid1,'PON');PON=netcdf.getVar(ncid1,varid,[0,0,0,t_start],[xd,yd,zd,num_days]);PON=double(PON);
varid=netcdf.inqVarID(ncid1,'OP');OP=netcdf.getVar(ncid1,varid,[0,0,0,t_start],[xd,yd,zd,num_days]);OP=double(OP);
varid=netcdf.inqVarID(ncid1,'sp');sp=netcdf.getVar(ncid1,varid,[0,0,0,t_start],[xd,yd,zd,num_days]);sp=double(sp);
varid=netcdf.inqVarID(ncid1,'lp');lp=netcdf.getVar(ncid1,varid,[0,0,0,t_start],[xd,yd,zd,num_days]);lp=double(lp);
varid=netcdf.inqVarID(ncid1,'sz');sz=netcdf.getVar(ncid1,varid,[0,0,0,t_start],[xd,yd,zd,num_days]);sz=double(sz);
varid=netcdf.inqVarID(ncid1,'lz');lz=netcdf.getVar(ncid1,varid,[0,0,0,t_start],[xd,yd,zd,num_days]);lz=double(lz);
varid=netcdf.inqVarID(ncid1,'pz');pz=netcdf.getVar(ncid1,varid,[0,0,0,t_start],[xd,yd,zd,num_days]);pz=double(pz);

% Get diagnostics
varid=netcdf.inqVarID(ncid2,'l_lim_sp');l_lim_sp=netcdf.getVar(ncid2,varid,[0,0,0,t_start-1],[xd,yd,zd,num_days]);l_lim_sp=double(l_lim_sp);
varid=netcdf.inqVarID(ncid2,'l_lim_lp');l_lim_lp=netcdf.getVar(ncid2,varid,[0,0,0,t_start-1],[xd,yd,zd,num_days]);l_lim_lp=double(l_lim_lp);
varid=netcdf.inqVarID(ncid2,'phy_chl');chl=netcdf.getVar(ncid2,varid,[0,0,0,t_start-1],[xd,yd,zd,num_days]);chl=double(chl);
varid=netcdf.inqVarID(ncid2,'THETA');tmp=netcdf.getVar(ncid2,varid,[0,0,0,t_start-1],[xd,yd,zd,num_days]);tmp=double(tmp);
netcdf.close(ncid1);netcdf.close(ncid2);netcdf.close(ncid3);

% Note: t_start-1 is used for diagnostics  because the model does not output an "initial condition" like it does for tracer fields
% Note: light limtation for small phytos and large phytos can be calculated outside the model but for convienence I have had the model output the light lim.
% Note: temperature is an input, but for this experiment I had it output temperature to make sure files were being read in correctly

% Rotate 3D fields  
NO2=[];NH2=[];SI2=[];DON2=[];PON2=[];OP2=[];sp2=[];lp2=[];
sz2=[];lz2=[];pz2=[];l_lim_sp2=[];l_lim_lp2=[];chl2=[];
l_lim_sp2=[];l_lim_lp2=[];chl2=[];tmp2=[];xc_all2=[];yc_all2=[];

xc=xc';
yc=yc';

for tt=1:size(NO,4)
for zz=1:size(NO,3)
NO2(:,:,zz,tt)=NO(:,:,zz,tt)';
NH2(:,:,zz,tt)=NH(:,:,zz,tt)';
SI2(:,:,zz,tt)=SI(:,:,zz,tt)';
DON2(:,:,zz,tt)=DON(:,:,zz,tt)';
PON2(:,:,zz,tt)=PON(:,:,zz,tt)';
OP2(:,:,zz,tt)=OP(:,:,zz,tt)';
sp2(:,:,zz,tt)=sp(:,:,zz,tt)';
lp2(:,:,zz,tt)=lp(:,:,zz,tt)';
sz2(:,:,zz,tt)=sz(:,:,zz,tt)';
lz2(:,:,zz,tt)=lz(:,:,zz,tt)';
pz2(:,:,zz,tt)=pz(:,:,zz,tt)';
l_lim_sp2(:,:,zz,tt)=l_lim_sp(:,:,zz,tt)';
l_lim_lp2(:,:,zz,tt)=l_lim_lp(:,:,zz,tt)';
chl2(:,:,zz,tt)=chl(:,:,zz,tt)';
tmp2(:,:,zz,tt)=tmp(:,:,zz,tt)';
end
end
NO=NO2;clear NO2;NH=NH2;clear NH2;SI=SI2;clear SI2;DON=DON2;clear DON2;
PON=PON2;clear PON2;OP=OP2;clear OP2;sp=sp2;clear sp2;
lp=lp2;clear lp2;sz=sz2;clear sz2;lz=lz2;clear lz2;pz=pz2;clear pz2;
l_lim_sp=l_lim_sp2;clear l_lim_sp2;l_lim_lp=l_lim_lp2;clear l_lim_lp2;
chl=chl2;clear chl2;
tmp=tmp2;clear tmp2; 

% Combine temporarily into one array 
if pp==1
t1=1;
t2=size(NO,4);
else
t1=t2+1;
t2=t2+td;
end

xc_all(:,:,dd)=xc;
yc_all(:,:,dd)=yc;
NO_all(:,:,:,t1:t2,dd)=NO;
NH_all(:,:,:,t1:t2,dd)=NH;
SI_all(:,:,:,t1:t2,dd)=SI;
DON_all(:,:,:,t1:t2,dd)=DON;
PON_all(:,:,:,t1:t2,dd)=PON;
OP_all(:,:,:,t1:t2,dd)=OP;
sp_all(:,:,:,t1:t2,dd)=sp;
lp_all(:,:,:,t1:t2,dd)=lp;
sz_all(:,:,:,t1:t2,dd)=sz;
lz_all(:,:,:,t1:t2,dd)=lz;
pz_all(:,:,:,t1:t2,dd)=pz;
l_lim_sp_all(:,:,:,t1:t2,dd)=l_lim_sp;
l_lim_lp_all(:,:,:,t1:t2,dd)=l_lim_lp;
chl_all(:,:,:,t1:t2,dd)=chl;
tmp_all(:,:,:,t1:t2,dd)=tmp;

end % for pp=1:length(ptr_files)
toc

end % for dd=1:length(dir_names)

display('')
display('Organize Tiles Into Single Arrays')
display('')
tic

% Get lon and lat range for entrie domain
dom_x1=min(reshape(xc_all,[1,numel(xc_all)]));
dom_x2=max(reshape(xc_all,[1,numel(xc_all)]));
dom_y1=min(reshape(yc_all,[1,numel(yc_all)]));
dom_y2=max(reshape(yc_all,[1,numel(yc_all)]));

% Get lon and lat range for each tile
x_range(:,1)=squeeze(xc_all(1,1,:));
x_range(:,2)=squeeze(xc_all(1,end,:));
y_range(:,1)=squeeze(yc_all(1,1,:));
y_range(:,2)=squeeze(yc_all(end,1,:));
lon(:,1)=sort(unique(x_range(:,1)));
lon(:,2)=sort(unique(x_range(:,2)));
lat(:,1)=sort(unique(y_range(:,1)));
lat(:,2)=sort(unique(y_range(:,2)));

% Re-Arrange tiles into 3D field
NO_all2=[];NH_all2=[];SI_all2=[];DON_all2=[];PON_all2=[];
OP_all2=[];sp_all2=[];lp_all2=[];sz_all2=[];lz_all2=[];
pz_all2=[];l_lim_sp_all2=[];l_lim_lp_all2=[];chl_all2=[];
l_lim_sp_all2=[];l_lim_lp_all2=[];chl_all2=[];xc_all2=[];yc_all2=[];

t_xdim=size(xc,2);
t_ydim=size(xc,1);

y1=1;
y2=t_ydim;
for yy=1:size(lat,1)
x1=1;
x2=t_xdim;
for xx=1:size(lon,1)
ind=find(x_range(:,1)==lon(xx,1) & x_range(:,2)==lon(xx,2) & y_range(:,1)==lat(yy,1) & y_range(:,2)==lat(yy,2));
xc_all2(y1:y2,x1:x2)=xc_all(:,:,ind);
yc_all2(y1:y2,x1:x2)=yc_all(:,:,ind);
NO_all2(y1:y2,x1:x2,:,:)=NO_all(:,:,:,:,ind);
NH_all2(y1:y2,x1:x2,:,:)=NH_all(:,:,:,:,ind);
SI_all2(y1:y2,x1:x2,:,:)=SI_all(:,:,:,:,ind);
DON_all2(y1:y2,x1:x2,:,:)=DON_all(:,:,:,:,ind);
PON_all2(y1:y2,x1:x2,:,:)=PON_all(:,:,:,:,ind);
OP_all2(y1:y2,x1:x2,:,:)=OP_all(:,:,:,:,ind);
sp_all2(y1:y2,x1:x2,:,:)=sp_all(:,:,:,:,ind);
lp_all2(y1:y2,x1:x2,:,:)=lp_all(:,:,:,:,ind);
sz_all2(y1:y2,x1:x2,:,:)=sz_all(:,:,:,:,ind);
lz_all2(y1:y2,x1:x2,:,:)=lz_all(:,:,:,:,ind);
pz_all2(y1:y2,x1:x2,:,:)=pz_all(:,:,:,:,ind);
l_lim_sp_all2(y1:y2,x1:x2,:,:)=l_lim_sp_all(:,:,:,:,ind);
l_lim_lp_all2(y1:y2,x1:x2,:,:)=l_lim_lp_all(:,:,:,:,ind);
chl_all2(y1:y2,x1:x2,:,:)=chl_all(:,:,:,:,ind);
tmp_all2(y1:y2,x1:x2,:,:)=tmp_all(:,:,:,:,ind);

x1=x2+1;
x2=x2+t_xdim;
end
y1=y2+1;
y2=y2+t_ydim;
end
NO_all=NO_all2;clear NO_all2;NH_all=NH_all2;clear NH_all2;
SI_all=SI_all2;clear SI_all2;DON_all=DON_all2;clear DON_all2;
PON_all=PON_all2;clear PON_all2;OP_all=OP_all2;clear OP_all2;
sp_all=sp_all2;clear sp_all2;lp_all=lp_all2;clear lp_all2;
sz_all=sz_all2;clear sz_all2;lz_all=lz_all2;clear lz_all2;
pz_all=pz_all2;clear pz_all2;
l_lim_sp_all=l_lim_sp_all2;clear l_lim_sp_all2;xc_all=xc_all2;clear xc_all2;
l_lim_lp_all=l_lim_lp_all2;clear l_lim_lp_all2;yc_all=yc_all2;clear yc)all2;
chl_all=chl_all2;clear chl_all2;tmp_all=tmp_all2;clear tmp_all2;

toc
display('Done Reading State Variables')

% Set Land Vals
NO_all(find(NO_all==0))=NaN;
NH_all(find(NH_all==0))=NaN;
SI_all(find(SI_all==0))=NaN;
DON_all(find(DON_all==0))=NaN;
PON_all(find(PON_all==0))=NaN;
OP_all(find(OP_all==0))=NaN;
sp_all(find(sp_all==0))=NaN;
lp_all(find(lp_all==0))=NaN;
sz_all(find(sz_all==0))=NaN;
lz_all(find(lz_all==0))=NaN;
pz_all(find(pz_all==0))=NaN;
l_lim_sp_all(find(l_lim_sp_all==0))=NaN;
l_lim_lp_all(find(l_lim_lp_all==0))=NaN;
tmp_all(find(tmp_all==0))=NaN;
chl_all(find(chl_all==0))=NaN;

% Example Diagnostics and Figures
NO_field=mean(NO_all,4);
chl_field=mean(chl_all,4);

figure(1)
pcolor(NO_field(:,:,1));shading flat;colorbar;caxis([0 1])
title('Nitrate (mmol N m^3)')

figure(2)
pcolor(chl_field(:,:,1));shading flat;colorbar;caxis([0 1])
title('Chl (mg chl m^3)')











