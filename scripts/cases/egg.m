mrstModule add ad-blackoil ad-core mrst-gui ad-props deckformat blackoil-sequential
mrstModule add linearsolvers agmg blackoil-sequential incomp
mrstVerbose true
%% egg


[G, rock, fluid, deck] = setupEGG('realization', 0);
G = computeBoundingBoxes(G);
%rock.perm = mean(rock.perm).*ones(G.cells.num, 3);
[state0, model, schedule, nonlinear] = initEclipseProblemAD(deck, 'G', G, 'TimestepStrategy', 'none');

%% refine

%% only valid for egg model
gridfromdeck=true;
refine=[4,4,2];
deck_new = rmfield(deck,'SCHEDULE');
physdims=deck_new.GRID.cartDims.*[deck_new.GRID.DX(1),deck_new.GRID.DY(1),deck_new.GRID.DZ(1)]
if(gridfromdeck)
    grdecl = simpleGrdecl(deck_new.GRID.cartDims, 0, 'flat',true','undisturbed',true,'physdims',physdims);
    grdecl.ZCORN=min(deck.GRID.TOPS)+grdecl.ZCORN;
    deck_new.GRID=rmfield(deck_new.GRID,{'DX','DY','DZ','TOPS'})
    deck_new.GRID.ZCORN=grdecl.ZCORN;
    deck_new.GRID.COORD=grdecl.COORD;
    deck_new = refineDeck(deck_new,refine)
    %%
    G_new = initEclipseGrid(deck_new);
    G_new = computeGeometry(G_new);
    rock_new  = initEclipseRock(deck_new);
    rock_new  = compressRock(rock_new, G_new.cells.indexMap);
    G_new=processGRDECL(deck_new.GRID);
else
    %% to do refine an mrst grid
    error('refinement of mrst grid not implemented');
    [G_new, rock_new] = deal(G, rock); % or something
    G_new=cartGrid(refine.*G.cartDims,physdims);
    G_new.nodes.coords(:,3)=G_new.nodes.coords(:,3)+min(deck.GRID.TOPS);
end
figure(2),clf,plotGrid(G_new)   
    %%
G_new=computeGeometry(G_new);
G_new = computeBoundingBoxes(G_new);
[schedule_new,maxperf] = makeNewSchedule(schedule,model.G, G_new,rock_new)

%%
nstep=10;
mkdir('tmp')
case_name='EGG'
outputprefix=fullfile(pwd(),'tmp',case_name);
mkdir(outputprefix)
model_new=model;
model_new.G=G_new;
if(nstep>0)
    schedule_new.step.control=schedule_new.step.control(1:nstep)
    schedule_new.step.val=schedule_new.step.val(1:nstep)
end
deck_new = model2Deck(model_new, schedule_new, 'deck', deck_new,'gridfromdeck',true)
deck_new.RUNSPEC.WELLDIMS(2)=100;
writeDeck(deck_new, outputprefix)
%
pth = getDatasetPath('egg');
%deckfile_org=fullfile(pth,'MRST','Egg_Model_ECL.DATA');
%deck_org=readEclipseDeck(deckfile_org);
mkdir('tmp_sims')
deckfile=fullfile(outputprefix,[case_name,'.DATA']);
%deck_new_sim = readEclipseDeck(deckfile);
%
outputdir=fullfile(pwd(),'tmp_sims',case_name);
simulator='/home/hnil/Documents/GITHUB/OPM/opm_source/master/builds/release_mpi/opm-simulators/bin/flow'
%%
%[wellsols, states, reports_opm, extra] =... 
runDeckOPM(deckfile,...
                                                    'outputdir',outputdir,...
                                                    'simulator',simulator,...
                                                    'force_timestep',false,...
                                                    'verbose',true,...
                                                    'no_output',false,...
                                                    'np',20,...
                                                    'lineartol',1e-3,...
                                                    'threads',1,...
                                                    'strongdefault',true);                                             
  plotWellSols(wellsols)
  %% run ion mrst
  run_mrst=true;
  if(run_mrst)
      tic
      lastn = maxNumCompThreads(2);
      %buildLinearSolvers(true)
      %[CXXFLAGS, LINK, LIBS] = setupMexOperatorBuildFlags()
      %buildMexOperators(true)
      useMex=true;
      deck_mrst=readEclipseDeck(deckfile);
      deck_mrst=convertDeckUnits(deck_mrst);
      %[state0_mrst, model_mrst, schedule_mrst, nonlinear_mrst] = initEclipseProblemAD(deck_mrst);
      [state0_mrst, model_mrst, schedule_mrst, nonlinear_mrst] = initEclipseProblemAD(deck_mrst,'TimestepStrategy', 'ds', 'useCPR', true, ...
          'useMex', useMex);
      %[wellsols, states, reports] = simulateScheduleAD(state0_mrst, model_mrst, schedule_mrst)
      %
      problem = packSimulationProblem(state0_mrst, model_mrst, schedule_mrst, 'egg', ...
          'Name', 'egg_test', 'nonlinearsolver', nonlinear_mrst);
      clearPackedSimulatorOutput(problem,'Prompt',false)
      simulatePackedProblem(problem, 'continueOnError', false);
      toc
  end
  
  %%
  [wellsols_mrst, states_mrst, reports_mrst] = getPackedSimulatorOutput(problem);
  plotWellSols(wellsols_mrst)
  
  