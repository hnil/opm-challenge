mrstModule add ad-blackoil ad-core mrst-gui ad-props deckformat blackoil-sequential
mrstModule add linearsolvers agmg blackoil-sequential incomp
mrstVerbose true
%% egg
opm = mrstPath('opm-tests');
assert(~isempty(opm), 'You must register https://github.com/opm/opm-tests as a module!');
run_mrst=true;
org_case_name='SPE9_CP'
np=1;%number of mpi in opm number of treads for mrst
for ref=1:1
refine=[1,1,1]*ref;
[deck, output] = getDeckOPMData('spe9', org_case_name);
case_name=[org_case_name,'_RX_',num2str(refine(1)),'_RY_',num2str(refine(2)),'_RZ_',num2str(refine(3))];
%%
G = initEclipseGrid(deck,'SplitDisconnected',false);
G = computeGeometry(G);    
rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);
[state0, model, schedule, nonlinear] = initEclipseProblemAD(deck, 'G', G, 'TimestepStrategy', 'none');
%% only valid for egg model
gridfromdeck=true;
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

[schedule_new,maxperf] = makeNewSchedule(schedule,model.G, G_new,rock_new);
disp('new wells calculated')
%%
nstep=-1;
mkdir('tmp')
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
end
%%
deckfile=fullfile(outputprefix,[case_name,'.DATA']);
outputdir=fullfile(pwd(),'tmp_sims',case_name);
simulator='/home/hnil/Documents/GITHUB/OPM/opm_source/master/builds/release_mpi/opm-simulators/bin/flow'
%deckfile='/data/hnil/GITHUB/opm-data/spe9/SPE9_CP.DATA'
[wellsols, states, reports_opm, extra] = runDeckOPM(deckfile,...
                                                    'outputdir',outputdir,...
                                                    'simulator',simulator,...
                                                    'force_timestep',false,...
                                                    'verbose',false,...
                                                    'no_output',false,...
                                                    'np',np,...
                                                    'openmp',1,...
                                                    'lineartol',1e-2,...
                                                    'strongdefault',true);                                             
  plotWellSols(wellsols)
  %% run ion mrst
  
  if(run_mrst)
      tic
      lastn = maxNumCompThreads(np)
      %buildLinearSolvers(true)
      %[CXXFLAGS, LINK, LIBS] = setupMexOperatorBuildFlags()
      %buildMexOperators(true)
      useMex=true;
      deck_mrst=readEclipseDeck(deckfile);
      deck_mrst=convertDeckUnits(deck_mrst);

      [state0_mrst, model_mrst, schedule_mrst, nonlinear_mrst] = initEclipseProblemAD(deck_mrst,'TimestepStrategy', 'ds', 'useCPR', true, ...
          'useMex', useMex);
      problem = packSimulationProblem(state0_mrst, model_mrst, schedule_mrst, 'egg', ...
          'Name', 'egg_test', 'nonlinearsolver', nonlinear_mrst);
      clearPackedSimulatorOutput(problem,'Prompt',false)
      simulatePackedProblem(problem, 'continueOnError', false);
      toc
      %%
    [wellsols_mrst, states_mrst, reports_mrst] = getPackedSimulatorOutput(problem);
    T_mrst = cumsum(schedule_mrst.step.val);
    T_opm=reports_opm.ReservoirTime;
    ws={wellsols_mrst, wellsols};
    ws=reMapWellsols(ws);
    times={T_mrst,T_opm};
    names={'mrst','opm'}
    plotWellSols(ws,times,'datasetnames', names)
  end
  
  %%
  rt=getReportTimings(reports_mrst)
  disp('******************************* MRST **************************')
  ff=fields(rt)
  for i=1:numel(ff)
    fprintf('%s \t %g \n', ff{i},sum([rt.(ff{i})]))
  end

  disp('*************************** OPM *******************************')
  command = ['tail -n 25 ', fullfile(outputdir,case_name),'.DBG', ' | head -n 16'];
  system(command);
  
  
  
  