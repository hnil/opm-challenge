mrstModule add ad-core ad-blackoil ad-props mrst-gui project-nccs co2lab project-multiresolution-ve matlab_bgl coarsegrid deckformat
physDims = [5000 5000 900]
res = [20 20 10]*2
G=cartGrid(res,physDims)
%% smaller model
grdecl=simpleGrdecl(G.cartDims,0,'flat',true,'physDims',physDims,'undisturbed',true);
perm = 100*ones(G.cells.num,1);

grdecl.PERMX=perm;
grdecl.PERMY=perm;
grdecl.PERMZ=perm;
grdecl.PORO=0.1*ones(G.cells.num,1);

coarsen = [5 5 5]
% use_org=prod(coarsen)==1;
% if(~use_org)
% grdecl_c=cutGrdecl(grdecl,[1 400;1 180;1 70]);
% %grdecl_c=cutGrdecl(grdecl,[1 40;1 18;1 7]);
% %%
% grdecl_c = coarseGrdecl(grdecl_c,coarsen,'only_grid',false)
% %%
% else
%    grdecl_c=grdecl; 
% end
grdecl_c  =grdecl;
%%
mrstModule add libgeometry
G_c = mprocessGRDECL(grdecl_c)
G_c = mcomputeGeometry(G_c);
rock_c        = grdecl2Rock(grdecl_c, G_c.cells.indexMap);
rock_c.perm   = convertFrom(rock_c.perm, darcy);
%clf,plotCellData(G_c,log10(rock.perm(:,3)/darcy)),view(3)
perm = rock_c.perm(:,3)/darcy;
clf,plotCellData(G_c,perm,perm>0.2),view(3)
plotCellData(G_c,perm,'FaceAlpha',0.1),view(3)
colorbar
%% net to gross
%grdecl.NTG = gdata.NG;

%% make deck
% use SPE1CASE2 2P model fluid.
%deck_org = readEclipseDeck("//data/hnil/GITHUB/OPM/opm-tests/spe1/SPE1CASE2_2P.DATA");
deck_org = readEclipseDeck("/data/hnil/GITHUB/OPM/opm-tests/aquifer-oilwater/2D_OW_CTAQUIFER.DATA")
deck_org.SOLUTION=rmfield(deck_org.SOLUTION,'AQUCT');
deck_org.SOLUTION=rmfield(deck_org.SOLUTION,'AQUANCON');
deck_org.SOLUTION.EQUIL=[0 270 700 0 0 0 0 0 0 1 0]
%
deck=deck_org;
deck.GRID = grdecl_c;
%
deck.RUNSPEC.cartDims = deck.GRID.cartDims;
deck.RUNSPEC.DIMENS = deck.GRID.cartDims;
%
% set up wells

%
%x = linespace(
%path = 
wgen={}
path = [3500 400 400; 4000 400 400] 
wgens{1}=struct('wname','P1','sign',-1,'type','bhp','val',200*barsa,'refDepth',0,'path',path,'comp_i',[0 1 0]);
path = [3700 1000 400; 3700 1000 600] 
wgens{2}=struct('wname','I1','sign',1,'type','bhp','val',400*barsa,'refDepth',0,'path',path,'comp_i',[1 0 0]);
Wo=[];Ws=[]
for i=1:numel(wgens)
    wgen=wgens{i};
    swellpath = makeSingleWellpath(wgen.path);
    wellpath = combineWellPaths({swellpath});   
    Wo  = getWellFromPath(Wo, G_c, rock_c, wellpath, ...
                'comp_i', wgen.comp_i,...
                'val',wgen.val,...
                'type', wgen.type,...
                'sign', wgen.sign,...
                'Name', wgen.wname,...
                'refDepth',wgen.refDepth);        
%
    addpath('/data/hnil/BITBUCKET/mrst-bitbucket/projects/project-nccs/nccs2019/matlab/tofrancesca/')
    tt = computeTraversedCells(G_c, wgen.path)
    V = tt.coord2 - tt.coord1;
    Ws = addWell(Ws,G_c,rock_c,tt.cells,'lineSegments', V,...
                    'comp_i', wgen.comp_i,...
                    'val',wgen.val,...
                    'type', wgen.type,...
                    'sign', wgen.sign,...
                    'Name', wgen.wname,...
                    'refDepth',wgen.refDepth);   
end
%%
max_perf=0;
for i=1:numel(W)
    max_perf=max(max_perf,numel(W(i).cells))
end
%%

%
perm = rock_c.perm(:,3)/darcy;
clf,plotCellData(G_c,perm,perm>0.2,'FaceAlpha',0.5),view(3)
plotCellData(G_c,perm,'FaceAlpha',0.1),view(3)
plotWell(G_c,Wo)
colorbar
xlabel('x')
ylabel('y')
%
W=Ws;
deck.SCHEDULE.control(1)=mrstWellToControl(W, G_c, deck.SCHEDULE.control(1),deck.RUNSPEC,'add_wellindex',true)
steps=linspace(0,10,40)*365;
steps=diff(steps);
deck.SCHEDULE.step.control = ones(size(steps));
deck.SCHEDULE.step.val = steps;
deck.RUNSPEC.WELLDIMS(2)=max_perf;
deck.RUNSPEC.UNIFOUT=1;
deck.GRID.INIT=1;
deck.SCHEDULE.RPTSCHED=1;

run_mrst=false

%%
if(~use_org)
    case_name = ['CARTGRID','_',num2str(coarsen(1)),'_',num2str(coarsen(2)),'_',num2str(coarsen(3))];
else
    case_name = ['CARTGRID_ORG']
end
nw=1;
for i=1:numel(ws)
    nw = max([nw,numel(ws)]);
end
deck.RUNSPEC.WELLDIMS(2)=nw
outputprefix=fullfile(pwd(),'tmp_data',case_name);
tic;writeDeck(deck,outputprefix);toc;
deckfile=fullfile(outputprefix,[case_name,'.DATA']);
disp('finnished writing deck file');
deck_mrst_org = readEclipseDeck(deckfile);
deck_mrst = convertDeckUnits(deck_mrst_org)
%%
wellSols= {};states={};report={};
if(run_mrst)  
    [state0_mrst, model_mrst, schedule_mrst, nonlinear] = initEclipseProblemAD(deck,'useMexGeometry', false, 'useMexProcessGrid', false)       %
    [wellSols{1}, states{1}, report{1}]  = simulateScheduleAD(state0_mrst, model_mrst, schedule_mrst)
    %%
    [state0_mrst, model_mrst, schedule_mrst, nonlinear] = initEclipseProblemAD(deck_mrst_org,'useMexGeometry', false, 'useMexProcessGrid', false)       %
    [wellSols{2}, states{2}, report{2}]  = simulateScheduleAD(state0_mrst, model_mrst, schedule_mrst)
    %%    
    %plotToolbar(model_mrst.G,states{11)
    %plotWell(G_c,Wo)
    plotWellSols(wellSols)
end
%%
outputdir=fullfile(pwd(),'tmp_sims',case_name);
simulator='/home/hnil/Documents/GITHUB/OPM/opm_source/master/builds/release_mpi/opm-simulators/bin/flow'
[wellsols, states, reports_opm, extra] = runDeckOPM(deckfile,'outputdir',outputdir,'simulator',simulator,'force_timestep',false,'verbose',false,'no_output',false)
%%
if(run_mrst)
    wellSols{3}=wellsols;
    states{3}=states;
    all_wellsols=wellSols;
    [all_wellsols,ind]=reMapWellsols(all_wellsols);
    all_reports={report{1},report{2},reports_opm};
    times=cellfun(@(x) reshape(x.ReservoirTime,[],1), all_reports,'UniformOut',false);
    plotWellSols(all_wellsols,times,'datasetnames',{'mrst1','mrst2','opm'})
end




