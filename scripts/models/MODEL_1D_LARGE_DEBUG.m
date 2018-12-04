function [deck,case_name] =MODEL_1D_LARGE_DEBUG(moduledir,tfac)
mydir=fullfile(pwd,'data');
case_name='model_1d_large_debug';
if(tfac~=1)
    case_name = [upper(case_name),'num2str(tfac)'];
else
    case_name = [upper(case_name)];
end
gravity on 
Lx=2000;
H = 20;
grdecl=simpleGrdecl([50, 1, 1],0,'flat',true,'undisturbed', false)
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
%% assume template have empty regions
W=[];
wname = 'Inj';
sign = 1;
rateval = sum(g.cells.volumes*0.1)/(10*year);
type = 'rate';val = rateval;
%type = 'bhp';val = 300*barsa;
dd=min(g.cells.centroids(:,3));
%W=addWell(W,g,rock,1,'comp_i', [1 0 0], 'val', val, 'type', type, 'sign', sign, 'Name', wname,'refDepth',dd);
W=verticalWell(W,g,rock,1,1,1:g.cartDims(3),'comp_i', [1 0 0], 'val', val, 'type', type, 'sign', sign, 'Name', wname,'refDepth',dd);

%W(end).status=0;
wname = 'Prod';
sign = -1;
type = 'bhp';val = 180*barsa;
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
type = 'bhp';val = 370*barsa;

dd=min(g.cells.centroids(:,3));
W2=addWell(W2,g,rock,1,'comp_i', [1 0 0], 'val', val, 'type', type, 'sign', sign, 'Name', wname,'refDepth',dd);
%W(end).status=0;
wname = 'Prod';
sign = -1;
type = 'bhp';val = 150*barsa;
%type = 'orat';val = -0.004
%type = 'resv';val = 

W2=addWell(W2,g,rock,g.cells.num,'comp_i', [0 1 0], 'val', val, 'type', type, 'sign', sign, 'Name', wname,'refDepth',dd);
for i=1:numel(W)
 W2(i).status=true;
 %W(i).lims=[];
end

W3=[];
wname = 'Inj';
sign = 1;
%type = 'rate';val = sum(g.cells.volumes*0.1)/(10*year);
type = 'bhp';val = 370*barsa;

dd=min(g.cells.centroids(:,3));
W3=addWell(W3,g,rock,1,'comp_i', [1 0 0], 'val', val, 'type', type, 'sign', sign, 'Name', wname,'refDepth',dd);
%W(end).status=0;
wname = 'Prod';
sign = -1;
%type = 'bhp';val = 150*barsa;
type = 'orat';val = -0.002
%type = 'resv';val = 

W3=addWell(W3,g,rock,g.cells.num,'comp_i', [0 1 0], 'val', val, 'type', type, 'sign', sign, 'Name', wname,'refDepth',dd);
for i=1:numel(W)
 W3(i).status=true;
 %W(i).lims=[];
end



%W(1).status=false;
%W(end).status=0;
%% the deck is used to define the fluid
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
ns=2
for i = 0:ns
    deck.SCHEDULE.control(i*2+1)=mrstWellToControl(W,g, deck.SCHEDULE.control(1),deck.RUNSPEC,'add_wellindex',true)
    deck.SCHEDULE.control(i*2+2)=mrstWellToControl(W2,g, deck.SCHEDULE.control(1),deck.RUNSPEC,'add_wellindex',true)
end
% set on til orat
%deck.SCHEDULE.control(i*2+2)=mrstWellToControl(W3,g, deck.SCHEDULE.control(1),deck.RUNSPEC,'add_wellindex',true)
nc=numel(deck.SCHEDULE.control);
%% set timesteps
deck.SCHEDULE.step.val=(diff(linspace(0,2,nc+1)*365)/10)*tfac;
deck.SCHEDULE.step.control=ones(size(deck.SCHEDULE.step.val));
deck.SCHEDULE.step.control=[1:nc];
deck.SCHEDULE.control=deck.SCHEDULE.control(1:max(deck.SCHEDULE.step.control));
wellsols={};states={};reports={};
%% mrst case
deck_mrst= convertDeckUnits(deck)
G= initEclipseGrid(deck_mrst);
G = computeGeometry(G);
rock = initEclipseRock(deck_mrst);
rock = compressRock(rock,G.cells.indexMap);
fluid = initDeckADIFluid(deck_mrst);
model = selectModelFromDeck(G, rock, fluid, deck_mrst);
schedule_mrst = convertDeckScheduleToMRST(model, deck_mrst);
regions = getInitializationRegionsDeck(model, deck_mrst);
[state0, pressures] = initStateBlackOilAD(model, regions);
state0.rs=fluid.rsSat(state0.pressure);
state0.rv=fluid.rvSat(state0.pressure);
%state0.rs=ones(3,1)*60;
deck.SOLUTION=mrstStateToSolution(G,state0);