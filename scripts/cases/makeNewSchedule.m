function  [schedule_new,maxperf] = makeNewSchedule(schedule,G,G_new,rock_new)
maxperf=0;
schedule_new=schedule;
for i=1:numel(schedule.control)
    W_new = [];
    W = schedule.control(i).W;
    W = addTrajectories(W, G, 2);
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
            wnew.compi=W(k).compi;
            W_new=[W_new;wnew];
        end
       % hack for not supported wells in writing
%        if(strcmp(W_new(end).type,'default') || strcmp(W_new(end).type,'rate'))
%            W_new(end).type='bhp';
%            W_new(end).status=0;
%        end
%        if(strcmp(lower(W_new(end).type),'resv_history'))
%            W_new(end).type='resv';
%         end
    end
    schedule_new.control(i).W=W_new;
    W_prev=W;
    W_prev_new=W_new;
end