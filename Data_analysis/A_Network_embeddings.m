clear all
close all
clc

% Use this script to generate figures similar to those in Figure 3C of:
%Luna, R., Li, J., Bauer, R., & van Leeuwen, C. (2023). Retinal waves in adaptive 
%rewiring networks orchestrate convergence and divergence in the visual system.
%Network Neuroscience



% This script generates figures of network embeddings, out-degree and
% in-degree histograms and adjacency matrices:

% Initial random network: Figures 1, 100, 10, 1000 and 2
% Figure 1: shows network embedding with out-degrees
% Figure 10: shows network embedding with in-degrees

% Evolved network: Figures 3, 300, 30, 3000 and 4
% Figure 3: shows network embedding with out-degrees
% Figure 30: shows network embedding with in-degrees



cd '..\scripts\Output'

fid = py.open('simulation_params.pckl','rb');
sim_params = py.pickle.load(fid);
tau = double(sim_params{1}{1});

% Select the following parameters as desired
k = 1; %network configuration (e.g. out of 150 possible ones)
pDist = 0.8; %probability for proximity-based rewiring
pBoost = 0.1; %calcium wave probability




% Which are wave initiator nodes?
fid = py.open('boosted_nodes.pckl','rb');
boosted_nodes = py.pickle.load(fid);
wave_initiator_node = boosted_nodes{k};
wave_initiator_node = double(wave_initiator_node{k})+1;



%% Load initial random network
fid = py.open('initials_nor.pckl','rb'); 
A = py.pickle.load(fid);
coords = A(1);
adj = A(2);


COORDS = double( coords{1} {k} );
xdata = COORDS(:,1)';
ydata = COORDS(:,2)';
ADJ_mat = double( adj{1} {0} ); %

G = digraph(ADJ_mat');


%out-degree
figure(1)
EdgeColors = normalize(G.Edges.Weight,'range');
EdgeColors = [EdgeColors,EdgeColors,EdgeColors];
p = plot(G,'EdgeColor',EdgeColors,'XData',xdata,'YData',ydata,'MarkerSize',5,'ArrowPosition',1);


wcc = centrality(G,'outdegree','Importance',G.Edges.Weight);

caxis manual
caxis([0 50]); 


colormap(cool)
cob = colorbar;
c.LineWidth = 3;

p.NodeCData = wcc;

set(gcf,'color','w');
axis square
box off
h = gca;
h.YAxis.Visible = 'off';
h.XAxis.Visible = 'off';




outdeg_hist = outdegree(G);
figure(100)
hist(outdeg_hist,[0 5 10 15 20 25 30 35 40 45 50 55 60 65 70],'BinCounts',[0 5 10 15 20 25 30 35 40 45 50 55 60 65 70])
h = findobj(gca,'Type','patch');
h.FaceColor = [0.5 0.5 0.5];
h.LineWidth = 1.5;

axis([0 70 0 60])

xlabel('Out-degree')
yticks([0:5:60])
yticklabels([0:5:60])
ylabel('Count')

set(gca,'linewidth',2)
set(gca,'FontSize',30)
set(gcf,'color','w');
box off

axis square





%in-degree
figure(10)
EdgeColors = normalize(G.Edges.Weight,'range');
EdgeColors = [EdgeColors,EdgeColors,EdgeColors];
p = plot(G,'EdgeColor',EdgeColors,'XData',xdata,'YData',ydata,'MarkerSize',5,'ArrowPosition',1);

wcc = centrality(G,'indegree','Importance',G.Edges.Weight);

caxis manual
caxis([0 50]); 


colormap(cool)
cob = colorbar;
c.LineWidth = 3;

p.NodeCData = wcc;

set(gcf,'color','w');
axis square
box off
h = gca;
h.YAxis.Visible = 'off';
h.XAxis.Visible = 'off';





indeg_hist = indegree(G);
figure(1000)
hist(indeg_hist,[0 5 10 15 20 25 30 35 40 45 50 55 60 65 70],'BinCounts',[0 5 10 15 20 25 30 35 40 45 50 55 60 65 70])
h = findobj(gca,'Type','patch');
h.FaceColor = [0.5 0.5 0.5];
h.LineWidth = 1.5;

axis([0 70 0 60])

xlabel('In-degree')
yticks([0:5:60])
yticklabels([0:5:60])
ylabel('Count')

set(gca,'linewidth',2)
set(gca,'FontSize',30)
set(gcf,'color','w');
box off

axis square







ax = figure(2)

b = imagesc(flipud(1-ADJ_mat));
set(b,'AlphaData',~isnan(ADJ_mat))

colormap(ax,gray)%


xticks([1:5:size(ADJ_mat,2)])
xticklabels([0:5:101])
xlabel('out-node')
    
yticks([1:5:size(ADJ_mat,1)])
yticklabels(fliplr([0:5:100]))
ylabel('in-node')

%title('Adjacency matrix') 

set(gca,'linewidth',2)
set(gca,'FontSize',20)
set(gcf,'color','w');


axis square





%% load evolved network

fid = py.open('boost_nor_MatPlusCoords.pckl','rb');  %  initials_nor.pckl
A = py.pickle.load(fid);
coords = A(1);
adj = A(2);



COORDS = double( coords{1} {k} );
xdata = COORDS(:,1)';
ydata = COORDS(:,2)';
ADJ_mat = double( adj{1} {tau}{0.5}{pDist}{pBoost}{k}{0} ); %A_matrices[tau][p_in][pDist][pBoost][k][r]

fid = py.open('boosted_nodes.pckl','rb');
boosted_nodes = py.pickle.load(fid);
pacemaker = boosted_nodes{k}

G = digraph(ADJ_mat');


%out-degree
figure(3)

EdgeColors = normalize(G.Edges.Weight,'range');
EdgeColors = [EdgeColors,EdgeColors,EdgeColors];
p = plot(G,'EdgeColor',EdgeColors,'XData',xdata,'YData',ydata,'MarkerSize',5,'ArrowPosition',1);


wcc = centrality(G,'outdegree','Importance',G.Edges.Weight);

caxis manual
caxis([0 50]);

colormap(cool)
cob = colorbar;
c.LineWidth = 3;

p.NodeCData = wcc;


set(gcf,'color','w');
axis square
box off
h = gca;
h.YAxis.Visible = 'off';
h.XAxis.Visible = 'off';




outdeg_hist = outdegree(G);
figure(300)
hist(outdeg_hist,[0 5 10 15 20 25 30 35 40 45 50 55 60 65 70],'BinCounts',[0 5 10 15 20 25 30 35 40 45 50 55 60 65 70])
h = findobj(gca,'Type','patch');
h.FaceColor = [0.5 0.5 0.5];
h.LineWidth = 1.5;

axis([0 70 0 60])

xlabel('Out-degree')
yticks([0:5:60])
yticklabels([0:5:60])
ylabel('Count')

set(gca,'linewidth',2)
set(gca,'FontSize',30)
set(gcf,'color','w');
box off

axis square





%in-degree
figure(30)

EdgeColors = normalize(G.Edges.Weight,'range');
EdgeColors = [EdgeColors,EdgeColors,EdgeColors];
p = plot(G,'EdgeColor',EdgeColors,'XData',xdata,'YData',ydata,'MarkerSize',5,'ArrowPosition',1);

wcc = centrality(G,'indegree','Importance',G.Edges.Weight);

caxis manual
caxis([0 50]);

colormap(cool)
cob = colorbar;
c.LineWidth = 3;

p.NodeCData = wcc;


set(gcf,'color','w');
axis square
box off
h = gca;
h.YAxis.Visible = 'off';
h.XAxis.Visible = 'off';





indeg_hist = indegree(G);
figure(3000)
hist(indeg_hist,[0 5 10 15 20 25 30 35 40 45 50 55 60 65 70],'BinCounts',[0 5 10 15 20 25 30 35 40 45 50 55 60 65 70])
h = findobj(gca,'Type','patch');
h.FaceColor = [0.5 0.5 0.5];
h.LineWidth = 1.5;

axis([0 70 0 60])

xlabel('In-degree')
yticks([0:5:60])
yticklabels([0:5:60])
ylabel('Count')

set(gca,'linewidth',2)
set(gca,'FontSize',30)
set(gcf,'color','w');
box off

axis square





ax = figure(4)

b = imagesc(flipud(1-ADJ_mat));
set(b,'AlphaData',~isnan(ADJ_mat))


colormap(ax,gray)

xticks([1:5:size(ADJ_mat,2)])
xticklabels([0:5:101])
xlabel('out-node')
    
yticks([1:5:size(ADJ_mat,1)])
yticklabels(fliplr([0:5:100]))
ylabel('in-node')


set(gca,'linewidth',2)
set(gca,'FontSize',20)
set(gcf,'color','w');


axis square



