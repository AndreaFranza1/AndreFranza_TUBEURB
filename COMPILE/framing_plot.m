function [choice] = framing_plot( Nodes_cordinates, Element_nodes )
%UNTITLED2 plots 2D frame with nodes and element numbering
%   Detailed explanation goes here

%plotting (it is prefered to be isulated as a seprate function)
%Adjacency matrix
adj=zeros(length(Nodes_cordinates),length(Nodes_cordinates));
for j=1:length(Element_nodes)
    adj(Element_nodes(j,1),Element_nodes(j,2))=1;
    adj(Element_nodes(j,2),Element_nodes(j,1))=1;
end

figure('name','Geometry');
hold on
gplot(adj,Nodes_cordinates)
% plot(Nodes_cordinates(:,1),Nodes_cordinates(:,2),'o');
set(gca, 'visible', 'off') ;

%Nodes numbering
for j = 1:length(Nodes_cordinates);
    text(Nodes_cordinates(j,1),Nodes_cordinates(j,2),int2str(j),'FontSize',8);
end

%member numbering
for j=1:length(Element_nodes);
    x05=mean(Nodes_cordinates(Element_nodes(j,:),1));
    y05=mean(Nodes_cordinates(Element_nodes(j,:),2));
    text(x05,y05,int2str(j),'FontSize',8,'color','r');
end

choice1 = questdlg('Geometry is OK?','Geometry visual check','Yes','No','Yes');
choice=strcmp(choice1,'No');
if choice==1
    display('Frame geometry need to be revised')
    text(1,-1,'Frame geometry need to be revised')
    return
else
%     display('Frame geometry is OK')
    close Geometry
end

