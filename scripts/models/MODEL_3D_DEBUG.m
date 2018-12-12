function [deck,case_name] =MODEL_1D_DEBUG(moduledir,layers)
case_name=['model_3d_debug_',num2str(layers)];
case_name = upper(case_name);
gravity on 
Lx=2000;
H = 20;
grdecl=simpleGrdecl([10, 10, layers],0,'flat',true,'undisturbed', false)
grdecl.ZCORN=grdecl.ZCORN*H*4+1000;
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
type = 'bhp';val = 320*barsa;
dd=min(g.cells.centroids(:,3));
%W=addWell(W,g,rock,1,'comp_i', [1 0 0], 'val', val, 'type', type, 'sign', sign, 'Name', wname,'refDepth',dd);
mid = floor(g.cartDims/2)+1;
W=verticalWell(W,g,rock,mid(1),mid(2),1,'comp_i', [1 0 0], 'val', val, 'type', type, 'sign', sign, 'Name', wname,'refDepth',dd);

%W(end).status=0;
wname = 'Prod1';
sign = -1;
type = 'bhp';val = 220*barsa;
%type = 'orat';val = -0.004
%type = 'resv';val = -0.004
W=verticalWell(W,g,rock,1,1,1:g.cartDims(3),'comp_i', [0 1 0], 'val', val, 'type', type, 'sign', sign, 'Name', wname,'refDepth',dd);


%%{
wname = 'Prod3';
W=verticalWell(W,g,rock,g.cartDims(1),g.cartDims(2),1:g.cartDims(3),'comp_i', [0 1 0], 'val', val, 'type', type, 'sign', sign, 'Name', wname,'refDepth',dd);
%W(end).status=0;
%%{
wname = 'Prod4';
sign = -1;
%type = 'orat';val = -0.001
type = 'bhp';val = 250*barsa;
W=verticalWell(W,g,rock,1,g.cartDims(2),1:g.cartDims(3),'comp_i', [0 1 0], 'val', val, 'type', type, 'sign', sign, 'Name', wname,'refDepth',dd);
wname = 'Prod2'
W=verticalWell(W,g,rock,g.cartDims(1),1,1:g.cartDims(3),'comp_i', [0 1 0], 'val', val, 'type', type, 'sign', sign, 'Name', wname,'refDepth',dd);
%}
for i=1:numel(W)
 W(i).status=true;
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
ns=0
for i = 0:ns
    deck.SCHEDULE.control(i*2+1)=mrstWellToControl(W,g, deck.SCHEDULE.control(1),deck.RUNSPEC,'add_wellindex',true)
    deck.SCHEDULE.control(i*2+2)=mrstWellToControl(W,g, deck.SCHEDULE.control(1),deck.RUNSPEC,'add_wellindex',true)
end
nc=numel(deck.SCHEDULE.control);
%% set timesteps
nc=numel(deck.SCHEDULE.control);
nc=1;
% set timesteps
n=15

%deck.SCHEDULE.step.val=[[1:4].^2*5,diff(linspace(0,2,n+1)*365)*10];
deck.SCHEDULE.step.val=[diff(linspace(0,2,n+1)*365)*5];
deck.SCHEDULE.step.control=ones(size(deck.SCHEDULE.step.val));
nsteps=numel(deck.SCHEDULE.step.control);
deck.SCHEDULE.step.control=[1:nc];
control=repmat([1:nc],[ceil(nsteps/nc),1]);
control=reshape(control,[],1);
deck.SCHEDULE.step.control=control(1:nsteps);
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
model.AutoDiffBackend = DiagonalAutoDiffBackend();
model.FacilityModel = UniformFacilityModel(model);
schedule_mrst = convertDeckScheduleToMRST(model, deck_mrst);
regions = getInitializationRegionsDeck(model, deck_mrst);
[state0, pressures] = initStateBlackOilAD(model, regions);
state0.rs=fluid.rsSat(state0.pressure);
state0.rv=fluid.rvSat(state0.pressure);
%state0.rs=ones(3,1)*60;
schedule_tmp =schedule_mrst;
schedule_tmp.step.val=[1:4].^2*2*day;
schedule_tmp.step.control=ones(size(schedule_tmp.step.val));
%%{
model.AutoDiffBackend = DiagonalAutoDiffBackend();
model.FacilityModel = UniformFacilityModel(model);
[wellsols_tmp, states_temp, reports_tmp] = simulateScheduleAD(state0, model, schedule_tmp);
state0=states_temp{end};
deck.SOLUTION=mrstStateToSolution(G,state0);
ns= numel(deck.SCHEDULE.step.control);
deck.SCHEDULE.control=repmat(deck.SCHEDULE.control(1),[ns,1])
deck.SCHEDULE.step.control=[1:ns]';
end