function stop = PlotFnForOptimiser(x,optimValues,state)

%note to self: need to have x as an input even if it's not used, otherwise the function won't run

persistent yFn it
stop = false;
switch state
    case 'iter'
        if optimValues.iteration == 0
            it = optimValues.iteration;            
                yFn = optimValues.fval;
                plot(it,yFn,'-','Marker','o');
                xlabel('Iteration'); ylabel('Objective function value')
        else
            it = [it optimValues.iteration];
                yFn = [yFn optimValues.fval];
                plot(it,yFn,'-','Marker','o');
                xlabel('Iteration'); ylabel('Objective function value')
        end
end

