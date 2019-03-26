

clear all;
close all;
addpath ~/Dropbox/0Work/MATLAB/m_map;


DD=load('oxygen_fk.txt');
[NOBS,dum] = size(DD);


xy=[-34,11,-32,-3];  
lat0=-17.5;

%domain[1]=xy[1]*111e3*cos(lat0*pi/180.0)
%domain[2]=xy[2]*111e3*cos(lat0*pi/180.0)
%domain[3]=xy[3]*111e3
%domain[4]=xy[4]*111e3


for i=1:NOBS
	DD(i,1) = DD(i,1)/(111e3 * cos(lat0*pi/180.0));	
	DD(i,2) = DD(i,2)/(111e3);
end





figure(1);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.5, 0.5, 0.5]);

m_proj('Mercator','lat',[-35 0],'lon',[-40 15]);
m_coast('patch',[.9 .9 .9],'edgecolor','none');
m_grid('linestyle','none','box','fancy','tickdir','out', 'fontsize', 20);
hold on;


p1=210;
p2=260;

cm=jet(2^8);
nc=size(cm,1);

colormap(cm);
cb = colorbar;
cb.FontSize = 20;
set(cb,'position',[0.86 .2 .02 .6])
caxis([p1,p2]);


for j=1:NOBS
	idx = ceil( (nc-1)* ( DD(j,3) - p1 )/(p2-p1) );
	m_plot(DD(j,1), DD(j,2), '.', 'markersize', 50, 'color', cm(idx,:));
end
m_plot([xy(1), xy(2), xy(2), xy(1), xy(1)], [xy(3), xy(3), xy(4), xy(4), xy(3)],'-k', 'linewidth', 4);
