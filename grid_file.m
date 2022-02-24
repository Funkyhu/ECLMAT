%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Generate grid file for ECLIPSE (structured grids) %%%%%
%%%%% Convert ECL grids to MODFLOW grids                %%%%%
%%%%% by Linwei Hu, @ Uni Kiel, 2017-08-22              %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc
clear
clf

%% Specify coordinate origin and discretization size 

ECL2MOD=0;    % 0: only ECLIPSE grid file is generated;  1: generate MODFLOW meshes also
ECL_GRID=0;   % 0: not generate ECLIPSE grid file;       1: ECLIPSE grid file

folder='G:\SNF\Eclipse\model_input\discretization\ECLIPSE\';
folder_check(folder)

folder1='G:\SNF\Eclipse\model_input\discretization\MODFLOW\';
folder_check(folder1)

% set properties of exported figures
figure_format='.emf';
resolution=1000;
printfig=0;    % print figure or not (1 or 0)

% specify coordinate origin
x0=0;      % left corner  (northwest)
y0=0;      % upper corner (northwest)
z0=1500;   % top elevation (in ECL it is positive, meaning z m below surface)
z_datum=z0;

dx=load(strcat(folder,'dx_test.txt'));
dy=load(strcat(folder,'dy_test.txt'));
dz=load(strcat(folder,'dz_test.txt'));

figfile1=strcat(folder,'grid_layer.emf'); delete_file(figfile1)
figfile2=strcat(folder,'ecl_vs_mod.emf'); delete_file(figfile2)
figfile3=strcat(folder,'grid.emf'); delete_file(figfile3)
gridfile=strcat(folder,'model_grid_test.inc'); delete_file(gridfile)


% specify z is negative or positive value
if z_datum>0
    z_index=1;     % -1: negative     1: positive
else
    z_index=-1;
end


%% Generate grid file for ECLIPSE
% corner points
xcoord=[x0]; ycoord=[y0]; zcoord=[z0];
for i=1:length(dx)
    x0=x0+dx(i);
    xcoord=[xcoord;x0];
end

for i=1:length(dy)
    y0=y0+dy(i);
    ycoord=[ycoord;y0];
end

for i=1:length(dz)
    z0=z0+dz(i);
    zcoord=[zcoord;z0];
end

NX=length(xcoord); nx=NX-1;
NY=length(ycoord); ny=NY-1;
NZ=length(zcoord); nz=NZ-1;

[XCOORD,YCOORD,ZCOORD]=meshgrid(xcoord,ycoord,zcoord);

save(strcat(folder,'x.txt'),'-ascii','xcoord');
save(strcat(folder,'y.txt'),'-ascii','ycoord');

% re-define ZCOORD (Gaussian distribution)
x_c=(xcoord(end)-xcoord(1))/2;
sigmax=300;
y_c=(ycoord(end)-ycoord(1))/2;
sigmay=300;
A=800;

[X,Y]=meshgrid(xcoord,ycoord);
if z_index==-1
    Z0=z_datum+A.*exp(-(X-x_c).^2./(2.*sigmax.^2)-(Y-y_c).^2./(2.*sigmay.^2));
elseif z_index==1
    Z0=z_datum-A.*exp(-(X-x_c).^2./(2.*sigmax.^2)-(Y-y_c).^2./(2.*sigmay.^2));
end
            % elevation at top layer
Z0(find(abs(Z0-z_datum)<0.001))=z_datum;

ZCOORD(:,:,1)=Z0; 

for i=1:nz
    ZCOORD(:,:,i+1)=ZCOORD(:,:,i)-dz(i);
end



% visualize the grid
gcf1=figure(1);
col_nr=5;
set(gcf1,'position',([5 50 2000 2000/(col_nr+1)]))

for i=1:NZ
    subplot(ceil(NZ/col_nr),col_nr,i)
    surf(XCOORD(:,:,i),YCOORD(:,:,i),ZCOORD(:,:,i))
    title(strcat('layer',32,num2str(i)))
    xlabel('X (m)');ylabel('Y (m)');zlabel('Z (m)');
    view(0,0)
    ZCOORD_btm=ZCOORD(:,:,NZ);
    ZCOORD_top=ZCOORD(:,:,1);
%     xlim([6000 10000])
    zlim([min(ZCOORD_btm(:)),max(ZCOORD_top(:))])
end


if printfig==1
    print_fig(figfile1,figure_format,resolution)
end


gcf2=figure(2);
for i=1:NZ
    h=mesh(X,Y,ZCOORD(:,:,i));hold on
    set(h,'edgecolor','k')
end
title('model discretization')
xlabel('X (m)');ylabel('Y (m)');zlabel('Z (m)');
% axis([0 16000 0 16000])
pbaspect([1 1 0.2])
view([0 -90])
hold off

if printfig==1
    print_fig(figfile3,figure_format,resolution)
end


% save the XYZ COORD data
XCOORD_total=[];
YCOORD_total=[];
ZCOORD_total=[];

for m=1:NZ
    XCOORD_total=[XCOORD_total;XCOORD(:,:,m)];
    YCOORD_total=[YCOORD_total;YCOORD(:,:,m)];
    ZCOORD_total=[ZCOORD_total;ZCOORD(:,:,m)];
end

% XYZ-coord (can be used for tecplot)
save(strcat(folder,'XCOORD_test.txt'),'-ascii','XCOORD_total');
save(strcat(folder,'YCOORD_test.txt'),'-ascii','YCOORD_total');
save(strcat(folder,'ZCOORD_test.txt'),'-ascii','ZCOORD_total');

% generate coordinate lines (COORD) and corner points (ZCORN)
Z1=reshape(ZCOORD(:,:,1),NX*NY,1);   
Z2=reshape(ZCOORD(:,:,NZ),NX*NY,1);   

X1=repmat(xcoord,NY,1);
Y1=repmat(ycoord,1,NX);
Y1=reshape(Y1',NX*NY,1);

COORD=[X1 Y1 Z1 X1 Y1 Z2];


ZCORN=[];
for k=1:nz
    % top
    for j=1:ny
        ZCOORD_k_T=ZCOORD(:,:,k);
        
        ZCOORD_k_T_j=ZCOORD_k_T(j,:);
        ZCOORD_k_T_j_1=ZCOORD_k_T(j+1,:);
        
        ZCOORD_k_T_j_NW=ZCOORD_k_T_j(1:NX-1);
        ZCOORD_k_T_j_NE=ZCOORD_k_T_j(2:NX);
        ZCOORD_k_T_j_SW=ZCOORD_k_T_j_1(1:NX-1);
        ZCOORD_k_T_j_SE=ZCOORD_k_T_j_1(2:NX);
        
        ZCORN_k_T_j=kron(ZCOORD_k_T_j_NW,[1 0])+kron(ZCOORD_k_T_j_NE,[0 1]);
        ZCORN_k_T_j_1=kron(ZCOORD_k_T_j_SW,[1 0])+kron(ZCOORD_k_T_j_SE,[0 1]);
        ZCORN=[ZCORN;ZCORN_k_T_j;ZCORN_k_T_j_1];
    end
    % bottom
    for j=1:ny
        ZCOORD_k_B=ZCOORD(:,:,k+1);
        
        ZCOORD_k_B_j=ZCOORD_k_B(j,:);
        ZCOORD_k_B_j_1=ZCOORD_k_B(j+1,:);
        
        ZCOORD_k_B_j_NW=ZCOORD_k_B_j(1:NX-1);
        ZCOORD_k_B_j_NE=ZCOORD_k_B_j(2:NX);
        ZCOORD_k_B_j_SW=ZCOORD_k_B_j_1(1:NX-1);
        ZCOORD_k_B_j_SE=ZCOORD_k_B_j_1(2:NX);
        ZCORN_k_B_j=kron(ZCOORD_k_B_j_NW,[1 0])+kron(ZCOORD_k_B_j_NE,[0 1]);
        ZCORN_k_B_j_1=kron(ZCOORD_k_B_j_SW,[1 0])+kron(ZCOORD_k_B_j_SE,[0 1]);
        ZCORN=[ZCORN;ZCORN_k_B_j;ZCORN_k_B_j_1];
    end
end


% 4. generate the grid file
if ECL_GRID==1
    
    fid=fopen(gridfile,'a+');
    
    fprintf(fid,'%s\n %s\n\n','PINCH',47);
    
    fprintf(fid,'%s\n','MAPUNITS');
    fprintf(fid,'   %s %s\n\n','METRES',47);
    
    fprintf(fid,'%s\n','GRIDUNIT');
    fprintf(fid,'   %s %s\n\n','METRES MAP',47);
    
    fprintf(fid,'%s\n','SPECGRID');
    fprintf(fid,'   %d %d %d %d %s %s\n\n',nx,ny,nz,1,'F',47);
    
    fprintf(fid,'%s\n','COORDSYS');
    fprintf(fid,'   %d %d %s\n\n',1,nz,47);
    
    fprintf(fid,'%s\n','COORD');
    fprintf(fid,'   %.2f %.2f %.2f %.2f %.2f %.2f\n',COORD');
    fprintf(fid,'   %s\n\n',47);
    
    fprintf(fid,'%s\n','ZCORN');
    fprintf(fid,'   %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n',ZCORN');
    fprintf(fid,'   %s\n\n',47);
    
    fprintf(fid,'%s\n','ACTNUM');
    fprintf(fid,'   %d%s%d\n\n',nx*ny*nz,'*',1);
    fprintf(fid,'   %s\n\n',47);
    
    fclose('all');
end


%%  ECLIPSE grid to MODELMUSE grid files %%
if ECL2MOD==1
    % corresponding MODFLOW model coordinates
    xcoord_c=zeros(1,length(xcoord)+1);     % xcoord for MODFLOW
    for i=1:length(xcoord)-1
        xcoord_c(i+1)=xcoord(i)+(xcoord(i+1)-xcoord(i))/2;
    end
    xcoord_c=xcoord_c';
    xcoord_c(1)=2*xcoord(1)-xcoord_c(2);
    xcoord_c(end)=2*xcoord(end)-xcoord_c(end-1);
    
    ycoord_c=zeros(1,length(ycoord)+1);     % ycoord for MODFLOW
    for i=1:length(ycoord)-1
        ycoord_c(i+1)=ycoord(i)+(ycoord(i+1)-ycoord(i))/2;
    end
    ycoord_c=ycoord_c';
    ycoord_c(1)=2*ycoord(1)-ycoord_c(2);
    ycoord_c(end)=2*ycoord(end)-ycoord_c(end-1);
    %     ycoord_c=ycoord_c*(-1);
    
    save(strcat(folder1,'x.txt'),'-ascii','xcoord_c');
    save(strcat(folder1,'y.txt'),'-ascii','ycoord_c');
    
    [XCOORDC,YCOORDC]=meshgrid(xcoord_c,ycoord_c);
    
    
    
    XCOORDC_total=[];
    YCOORDC_total=[];
    ZCOORDC_total=[];
    % zcoord for MODFLOW
    for i=1:NZ
        zz=ZCOORD(:,:,i);         % zcoord of ECLIPSE at each layer (top->btm, corner point)
        model_btm=ZCOORD(:,:,i);
        model_btm_c=interp2(XCOORD(:,:,i),YCOORD(:,:,i),model_btm,XCOORDC,YCOORDC,'linear'); % zcoord of MODFLOW at each layer (top->btm, center)
        model_BTM=model_btm_c;
        model_BTM(2:end,1:1)=zz(:,1);
        model_BTM(2:end,end:end)=zz(:,end);
        model_BTM(1:1,2:end)=zz(1,:);
        model_BTM(end:end,2:end)=zz(end,:);
        model_BTM(1,1)=(model_BTM(1,2)+model_BTM(2,1))/2;
        
        % compare with ECLIPSE
        gcf2=figure(2);
        set(gcf2,'position',([800 300 1000 400]))
        
        subplot(1,2,1)
        mesh(ZCOORD(:,:,i));view(0,0)
        hold on
        subplot(1,2,2)
        mesh(model_BTM);view(0,0)
        
        % save the XYZ COORD data (MODFLOW)
        
        
        XCOORDC_total=[XCOORDC_total;XCOORDC];
        YCOORDC_total=[YCOORDC_total;YCOORDC];
        ZCOORDC_total=[ZCOORDC_total;model_BTM];
        
        
        % XYZ-coord (can be used for tecplot)
        save(strcat(folder1,'XCOORD_test.txt'),'-ascii','XCOORDC_total');
        save(strcat(folder1,'YCOORD_test.txt'),'-ascii','YCOORDC_total');
        save(strcat(folder1,'ZCOORD_test.txt'),'-ascii','ZCOORDC_total');
        
        
        
        % save the data in MODELMUSE recognized format
        
        model_xcoord=reshape(XCOORDC,1,size(XCOORDC,1)*size(XCOORDC,2))';
        model_ycoord=reshape(YCOORDC,1,size(YCOORDC,1)*size(YCOORDC,2))';
        model_xy_btm=[model_xcoord model_ycoord (reshape(model_BTM,1,size(model_BTM,1)*size(model_BTM,2)))'];
        str2=strcat('model_XY_BTM',num2str(i-1),'=model_xy_btm;');eval(str2)   % XYZ coord for each layer (model_XY_BTM0 is the model_TOP)
        
        if i==1
            save(strcat(folder1,'Z_model_top.txt'),'-ascii','model_xy_btm');
        else
            save(strcat(folder1,'Z_model_btm',num2str(i-1),'.txt'),'-ascii','model_xy_btm');
        end
        
    end
    
    if printfig==1
        print_fig(figfile2,figure_format,resolution);
    end
    
end





fclose('all');
