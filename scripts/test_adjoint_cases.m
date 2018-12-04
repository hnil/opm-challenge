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

%mycases={'1D','1D_LARGE','3D'}
mycases={'1D'}
% make deck
for i=1:numel(mycases)
    mycase=mycases{i};
    switch mycase
        case '1D' 
            [deck,case_name] =MODEL_1D_DEBUG(moduledir)
        case '1D_LARGE'    
            [deck,case_name] =MODEL_1D_LARGE_DEBUG(moduledir,1)
        case '3D'    
            [deck,case_name] =MODEL_3D_DEBUG(moduledir,1)
        otherwise
            error()
    end
end
          
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


%[wellsols{1}, states{1}, reports{1}] = simulateScheduleAD(state0, model, schedule_mrst,'NonLinearSolver',nonlinear);
[wellsols{1}, states{1}, reports{1}] = simulateScheduleAD(state0, model, schedule_mrst);

npvopts   = {'OilPrice',             0.0*stb , ...
                 'GasPrice',             0.0 , ...
                 'GasInjectionCost',     0.0 , ...
                 'WaterProductionCost',  -1.0*stb() , ...
                 'WaterInjectionCost',   1.0*stb() , ...
                 'DiscountFactor',       0.0};                                  
wellsols_mrst=wellsols{1};
states_mrst = states{1};
%objh = @(tstep)NPVOW(G, wellsols_mrst, schedule_mrst, 'computePartials', true, 'tstep', tstep, npvopts{:});
obj = @(wellSols, schedule, varargin)  NPVBlackOil(model.G, wellSols, schedule, varargin{:}, npvopts{:});
objh = @(tstep)obj(wellsols_mrst, schedule_mrst, 'ComputePartials', true, 'tStep', tstep);
sensad = computeGradientAdjointAD(state0, states_mrst, model, schedule_mrst, objh);
                                     %'Parameters'    , params, ...
                                     %'ParameterTypes', paramTypes);
                                     

%%
sensad{:}
objad=obj(wellsols_mrst,schedule_mrst)
%{
obj_tmp =@(wellSols,states, schedule, varargin) obj(wellSols, schedule,varargin{:})
sensnum = computeGradientPerturbationAD(state0,  model, schedule_mrst, obj_tmp,'perturbation',1*barsa());
                                     %'Parameters'    , params, ...
                                     %'ParameterTypes', paramTypes);
sensnum{:}
schedule_mrst.step.val
%}
%%

%%
% run adjoint opm
tic
outputdir = fullfile(moduledir,'tmp_data',case_name);
simulator='/home/hnil/Documents/GITHUB/OPM/opm_source/devel_adjoint/builds/serial/opm-simulators/bin/flow'
[wellsols{2}, states{2}, reports{2}] = runDeckOPM(deckfile,'outputdir',outputdir,'simulator',simulator,'force_timestep',true,'do_adjoint',true,'verbose',true,'no_output',false);
copyfile(outputprefix,fullfile(outputdir,'inputfiles')) 
toc
opm_rep=reports{2};

wellsols=reMapWellsols(wellsols);
nstep = numel(schedule_mrst.step.val);
cderv = zeros(numel(schedule_mrst.control),numel(schedule_mrst.control(1).W))
ctrl =schedule_mrst.step.control;
for i=1:nstep
    cderv(ctrl(i),:)= cderv(ctrl(i),:) + opm_rep.adjoint.derv(i,:);
end

W= schedule_mrst.control(1).W;
[a,ind]=sort({W.name});
disp('********************************************************')
fprintf('objective\t mrst:\t');
fprintf('%f ',[objad{:}]);
fprintf('\n')
fprintf('objective\t obj:\t');
fprintf('%f ', opm_rep.adjoint.obj);
fprintf('\n');
for i=1:numel(sensad)
    fprintf('Well name \t \t')
    fprintf('%s\t\t',W.name)
    fprintf('\n')
    fprintf('Control %i :', i);
    fprintf('mrst der \t');
    fprintf('%f\t',sensad{i});
    fprintf('\n\t opm der \t');
    fprintf('%f\t', cderv(i,ind));
    fprintf('\n')
end
%times=cellfun(@(x) reshape(x.ReservoirTime,[],1), reports,'UniformOut',false);
%plotWellSols(wellsols, times,'datasetnames',{'mrst','flow test'},'field','qWs','wells','all')


%
%sensad{1}(2)
% prepear for opm run
%figure(33)
%clf,plotToolbar(G,states{1})
%view(3)
%
mrst_der=horzcat(sensad{:})';
opm_der=-cderv(:,ind);
max(abs(mrst_der(:)-opm_der(:))./opm_der(:))
[mrst_der-opm_der]
opm_der
%%
return
