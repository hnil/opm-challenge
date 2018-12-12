function [deck,case_name] =MODEL_1D_DEBUG(moduledir)
mydir=fullfile(pwd,'data');
case_name='model_1d_debug';
case_name = upper(case_name);
gravity on 
Lx=2000;
H = 20;
grdecl=simpleGrdecl([3, 1, 1],0,'flat',true,'undisturbed', false)
grdecl.ZCORN=grdecl.ZCORN*H+1000;
grdecl.COORD=grdecl.COORD*Lx;
g=processGRDECL(grdecl);
g=computeGeometry(g);
PERMX=100*ones(prod(g.cartDims),1);
PORO=0.1*ones(prod(g.cartDims),1);
perm = [PERMX(g.cells.indexMap),PERMX(g.cells.indexMap),0.1*PERMX(g.cells.indexMap)];

rock=struct('perm',perm,'poro',PORO(g.cells.indexMap));
rock.perm=rock.perm*darcy;
max_coord=max(g.nodes.coords,[],1);
min_coord=min(g.nodes.coords,[],1);
z_mean = mean(g.cells.centroids(:,3));
%plotGrid(g)
%grdecl_filename = fullfile(mydir,[case_name,'.GRDECL'])
%writeGRDECL(grdecl,grdecl_filename)
%%
% change runspec
%% keep props
% deck.PROPS
%% assume template have empty regions
W=[];
wname = 'Inj';
sign = 1;
type = 'rate';val = sum(g.cells.volumes*0.1)/(10*year);
%type = 'bhp';val = 370*barsa;

dd=min(g.cells.centroids(:,3));
W=addWell(W,g,rock,1,'comp_i', [1 0 0], 'val', val, 'type', type, 'sign', sign, 'Name', wname,'refDepth',dd);
%W(end).status=0;
wname = 'Prod';
sign = -1;
type = 'bhp';val = 150*barsa;
%type = 'orat';val = -0.004
%type = 'resv';val = -0.004

W=addWell(W,g,rock,g.cells.num,'comp_i', [0 1 0], 'val', val, 'type', type, 'sign', sign, 'Name', wname,'refDepth',dd);
for i=1:numel(W)
 W(i).status=true;
 %W(i).lims=[];
end

W2=[];
wname = 'Inj';
sign = 1;
%type = 'rate';val = sum(g.cells.volumes*0.1)/(10*year);
type = 'bhp';val = 350*barsa;

dd=min(g.cells.centroids(:,3));
%dd=0;
W2=addWell(W2,g,rock,1,'comp_i', [1 0 0], 'val', val, 'type', type, 'sign', sign, 'Name', wname,'refDepth',dd);
%W(end).status=0;
wname = 'Prod';
sign = -1;
ype = 'bhp';val = 153*barsa;
%type = 'orat';val = -0.004
%type = 'resv';val = -0.004

W2=addWell(W2,g,rock,g.cells.num,'comp_i', [0 1 0], 'val', val, 'type', type, 'sign', sign, 'Name', wname,'refDepth',dd);
for i=1:numel(W)
 W(i).status=true;
 %W(i).lims=[];
end

%W(1).status=false;
%W(end).status=0;
%% the deck is used to define the fluid
%temp= fullfile(moduledir,'template/1D_CASE_OLYMPUS.DATA');
temp= fullfile(moduledir,'template/1D_CASE.DATA');
deck_org = readEclipseDeck(temp);
deck_org.SCHEDULE.control= deck_org.SCHEDULE.control(1)
deck = deck_org;
deck.RUNSPEC.DIMENS=grdecl.cartDims;
nc=prod(grdecl.cartDims);
grdecl.PERMX=100*ones(nc,1);
grdecl.PERMY=grdecl.PERMX;
grdecl.PERMZ=0.1*grdecl.PERMX;
grdecl.PORO=0.1*ones(nc,1);
%% change grid section to grdecl
nc = prod(grdecl.cartDims);
deck.GRID = grdecl;
%% 
% set equil_num
deck.SOLUTION.EQUIL(1)=z_mean-20;
deck.SOLUTION.EQUIL(3)=z_mean+40;
deck.SOLUTION.EQUIL(5)=z_mean-20;
%% add cartesian indexes for the wells
%deck.SCHEDULE.control(1)= mrstWellToControl(W,g, deck.SCHEDULE.control(1))
ns=-1
for i = 0:ns
    deck.SCHEDULE.control(i*2+1)=mrstWellToControl(W,g, deck.SCHEDULE.control(1),deck.RUNSPEC,'add_wellindex',true)
    deck.SCHEDULE.control(i*2+2)=mrstWellToControl(W2,g, deck.SCHEDULE.control(1),deck.RUNSPEC,'add_wellindex',true)
end
if(ns<0)
    deck.SCHEDULE.control=mrstWellToControl(W2,g, deck.SCHEDULE.control(1),deck.RUNSPEC,'add_wellindex',true)
else
    deck.SCHEDULE.control(i*2+3)=mrstWellToControl(W,g, deck.SCHEDULE.control(1),deck.RUNSPEC,'add_wellindex',true)
end
nc=numel(deck.SCHEDULE.control);
% set timesteps
n=2
deck.SCHEDULE.step.val=diff(linspace(0,2,n+1)*365)/4;
deck.SCHEDULE.step.control=ones(size(deck.SCHEDULE.step.val));
nsteps=numel(deck.SCHEDULE.step.control);
deck.SCHEDULE.step.control=[1:nc];
control=repmat([1:nc],[ceil(nsteps/nc),1]);
control=reshape(control,[],1);
%
deck.SCHEDULE.step.control=control(1:nsteps);
deck.SCHEDULE.control=deck.SCHEDULE.control(1:max(deck.SCHEDULE.step.control));
wellsols={};states={};reports={};

%deck.PROPS.SWOF{1}(:,end)=0.0;
%% mrst case
%{
deck.RUNSPEC.UNIFOUT=1;
deck.GRID.INIT=1;
deck.SCHEDULE.RPTSCHED=1;
outputprefix=fullfile(moduledir,'data',case_name);
writeDeck(deck,outputprefix);
deckfile=fullfile(outputprefix,[case_name,'.DATA']);
deck_new=readEclipseDeck(deckfile);
%}
end

%state0.rs=ones(3,1)*60;
%deck.SOLUTION=mrstStateToSolution(G,state0);
%deck.SOLUTION=rmfield(deck.SOLUTION,'EQUIL')
% PRVD PDEW