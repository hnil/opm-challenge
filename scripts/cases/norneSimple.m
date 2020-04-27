mrstModule add ad-blackoil ad-core mrst-gui ad-props deckformat blackoil-sequential
mrstModule add linearsolvers agmg blackoil-sequential incomp
mrstVerbose true
%% egg
opm = mrstPath('opm-tests');
assert(~isempty(opm), 'You must register https://github.com/opm/opm-tests as a module!');
[deck, output] = getDeckOPMData('norne', 'NORNE_ATW2013');
%%
deck.RUNSPEC=rmfield(deck.RUNSPEC,{'EQLOPTS','ENDSCALE','SATOPTS','GRIDOPTS'})
deck.GRID=rmfield(deck.GRID,{'FAULTS','MULTFLT','FLUXNUM','PINCH'})
deck.PROPS=rmfield(deck.PROPS,{'SWATINIT','EHYSTR'});%%,'SCALECRS'});
G = initEclipseGrid(deck,'SplitDisconnected',false);
G = computeGeometry(G);    
rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);
[state0, model, schedule, nonlinear] = initEclipseProblemAD(deck, 'G', G, 'TimestepStrategy', 'none');





%% refine

%% only valid for egg model
gridfromdeck=true;
refine=[1,1,1];
deck_new = rmfield(deck,'SCHEDULE');
deck_new = refineDeck(deck_new,refine)
G_new = initEclipseGrid(deck_new,'SplitDisconnected',false);
G_new = computeGeometry(G_new);
rock_new  = initEclipseRock(deck_new);
rock_new  = compressRock(rock_new, G_new.cells.indexMap);
%G_new=processGRDECL(deck_new.GRID);
%%
figure(2),clf,plotGrid(G_new)   
    %%
G_new=computeGeometry(G_new);
G_new = computeBoundingBoxes(G_new);
%%
[schedule_new,maxperf] = makeNewSchedule(schedule,model.G, G_new,rock_new)
disp('new wells calculated')
%%
nstep=-1;
mkdir('tmp')
case_name='NORNE'
outputprefix=fullfile(pwd(),'tmp',case_name);
mkdir(outputprefix)
model_new=model;
model_new.G=G_new;
if(nstep>0)  
    schedule_new.step.control=schedule_new.step.control(1:nstep)
    schedule_new.step.val=schedule_new.step.val(1:nstep)
    mcont=max(schedule_new.step.control);
    schedule_new.control=schedule_new.control(1:mcont);
end
%%
deck_new = model2Deck(model_new, schedule_new, 'deck', deck_new,'gridfromdeck',true)
deck_new.RUNSPEC.WELLDIMS(2)=maxperf;
writeDeck(deck_new, outputprefix)
%%
deckfile=fullfile(outputprefix,[case_name,'.DATA']);
outputdir=fullfile(pwd(),'tmp_sims',case_name);
simulator='/home/hnil/Documents/GITHUB/OPM/opm_source/master/builds/release_mpi/opm-simulators/bin/flow'
[wellsols, states, reports_opm, extra] = runDeckOPM(deckfile,...
                                                    'outputdir',outputdir,...
                                                    'simulator',simulator,...
                                                    'force_timestep',false,...
                                                    'verbose',false,...
                                                    'no_output',false,...
                                                    'np',10,...
                                                    'openmp',1,...
                                                    'strongdefault',true);                                             
  plotWellSols(wellsols)
  %% run ion mrst
  run_mrst=true;
  if(run_mrst)
      tic
      lastn = maxNumCompThreads(1)
      %buildLinearSolvers(true)
      %[CXXFLAGS, LINK, LIBS] = setupMexOperatorBuildFlags()
      %buildMexOperators(true)
      %%
      useMex=true;
      deck_mrst=readEclipseDeck(deckfile);
      deck_mrst=convertDeckUnits(deck_mrst);
      %[state0_mrst, model_mrst, schedule_mrst, nonlinear_mrst] = initEclipseProblemAD(deck_mrst);
      [state0_mrst, model_mrst, schedule_mrst, nonlinear_mrst] = initEclipseProblemAD(deck_mrst,'TimestepStrategy', 'ds', 'useCPR', true, ...
          'useMex', useMex);
      %%
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
  
  