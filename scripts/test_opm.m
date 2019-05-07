%%mrst
%moduledir='/home/hnil/Documents/BITBUCKET/mrst-bitbucket/projects/project-ensablecase/';
addpath(fullfile(pwd,'models'));
moduledir='/data/hnil/BITBUCKET/mrst-bitbucket/projects/opm-challenge/'
mrstModule add ad-props deckformat mrst-gui ad-core ad-blackoil
mrstModule add run-datasets wellpaths deckformat optimization
mrstModule add wellpaths
if exist('OCTAVE_VERSION', 'builtin') ~= 0
    mrstModule add octave
    run set_reasonable_octave_defaults;
end
force_timestep=true;
%mycases={'1D','1D_LARGE','3D'}
mycases={'1D'}
% make deck
for i=1:numel(mycases)
    mycase=mycases{i};
    switch mycase
        case '1D' 
            [deck,case_name] =MODEL_1D_DEBUG(moduledir)
        case '1D_LARGE'    
            [deck,case_name] =MODEL_1D_LARGE_DEBUG(moduledir,0.2)
        case '2D'
            [deck,case_name] =MODEL_2D_DEBUG(moduledir,3)
        case '3D'    
            [deck,case_name] =MODEL_3D_DEBUG(moduledir,3)
        otherwise
            error()
    end
end
ind=1;
deck.SCHEDULE.control = deck.SCHEDULE.control(ind);
deck.SCHEDULE.step.control = deck.SCHEDULE.step.control(ind)
deck.SCHEDULE.step.val = deck.SCHEDULE.step.val(ind)
%write deck
deck.RUNSPEC.UNIFOUT=1;
deck.GRID.INIT=1;
deck.SCHEDULE.RPTSCHED=1;
outputprefix=fullfile(moduledir,'tmp_models',case_name);
writeDeck(deck,outputprefix);
deckfile=fullfile(outputprefix,[case_name,'.DATA']);
disp('finnished writing deck file');


deck_mrst = readEclipseDeck(deckfile);
deck_mrst= convertDeckUnits(deck_mrst)
[state0, model, schedule_mrst, nonlinear] = initEclipseProblemAD(deck_mrst,'useMexGeometry', false, 'useMexProcessGrid', false)
% seems like undersaturated is a bit difficult
if(isfield(model.fluid,'rsSat'))
    state0.rs=model.fluid.rsSat(state0.pressure);
    state0.rv=model.fluid.rvSat(state0.pressure);    
else
    state0.rs =state0.pressure*0;
    state0.rv =state0.pressure*0;
end
deck.SOLUTION=mrstStateToSolution(model.G,state0);
writeDeck(deck,outputprefix);
%%
nonlinear.LinearSolver = selectLinearSolverAD(model,'BackslashThreshold',  1)
%
%[wellsols{1}, states{1}, reports{1}] = simulateScheduleAD(state0, model, schedule_mrst,'NonLinearSolver',nonlinear);
[wellsols{1}, states{1}, reports{1}] = simulateScheduleAD(state0, model, schedule_mrst,'NonLinearSolver', nonlinear)

                                     

%%
% run adjoint opm
tic
outputdir = fullfile(moduledir,'tmp_data',case_name);
simulator='/home/hnil/Documents/GITHUB/OPM/opm_source/master/builds/debug/opm-simulators/bin/flow'
%simulator='/home/hnil/Documents/GITHUB/OPM/opm_source/devel_adjoint_back/builds/debug/opm-simulators/bin/flow'
[wellsols{2}, states{2}, reports{2}, extra] = runDeckOPM(deckfile,'outputdir',outputdir,'simulator',simulator,'force_timestep',force_timestep,'do_adjoint',false,'verbose',true,'no_output',false);
copyfile(outputprefix,fullfile(outputdir,'inputfiles')) 
toc
opm_rep=reports{2};

wellsols_org=wellsols;
[wellsols,ind]=reMapWellsols(wellsols);

%%
times=cellfun(@(x) reshape(x.ReservoirTime,[],1), reports,'UniformOut',false);
plotWellSols_old(wellsols,times,'datasetnames',{'mrst','flow test'},'field','qWs','wells','all')
