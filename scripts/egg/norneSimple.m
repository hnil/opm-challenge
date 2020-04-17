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

schedule_new=schedule;
maxperf=0;
for i=1:numel(schedule.control)
    W_new = [];
    W = schedule.control(i).W;
    W = addTrajectories(W, model.G, 2);
    for k = 1:numel(W)
        w = W(k);
        newperfs=true;
        if(i>1)
            w_prev=W_prev(k);
            if(all(w.cells == w_prev.cells) &&  all(w.WI == w_prev.WI))
                newperfs=false;
            end
        end
        if(newperfs)
            tmp = computeTraversedCellsNew(G_new, w.trajectory);
            W_new   = addWell(W_new, G_new, rock_new, tmp.cell, 'name', w.name, 'type', w.type, 'sign', ...
                w.sign, 'val', w.val, 'refDepth', w.refDepth, ...
                'compi', w.compi, 'lims', w.lims, ...
                'lineSegments', bsxfun(@times, tmp.vec, tmp.weight) );
            maxperf=max(maxperf,numel(W_new(end).cells));
        else
            wnew=W_prev_new(k);
            wnew.type = W(k).type;
            wnew.val = W(k).val;
            wnew.cstatus=W(k).status;
            wnew.cstatus=W(k).compi;
            W_new=[W_new;wnew];
        end
    end
    schedule_new.control(i).W=W_new;
    W_prev=W;
    W_prev_new=W_new;
end
disp('new wells calculated')
%
nstep=-5;
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
deck_new.RUNSPEC.WELLDIMS(2)=100;
writeDeck(deck_new, outputprefix)
%
deckfile=fullfile(outputprefix,[case_name,'.DATA']);
outputdir=fullfile(pwd(),'tmp_sims',case_name);
simulator='/home/hnil/Documents/GITHUB/OPM/opm_source/master/builds/release_mpi/opm-simulators/bin/flow'
[wellsols, states, reports_opm, extra] = runDeckOPM(deckfile,...
                                                    'outputdir',outputdir,...
                                                    'simulator',simulator,...
                                                    'force_timestep',false,...
                                                    'verbose',false,...
                                                    'no_output',false,...
                                                    'np',1,...
                                                    'openmp',1);                                             
  plotWellSols(wellsols)
  %% run ion mrst
  run_mrst=true;
  if(run_mrst)
      tic
      lastn = maxNumCompThreads(1)
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
  
  