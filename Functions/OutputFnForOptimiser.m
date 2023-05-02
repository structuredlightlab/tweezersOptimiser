function stop = OutputFnForOptimiser(x,optimValues,state)

%This function generates output arrays for an optimiser to track its progress.
%INPUT
% x                     values of the optimisation variables at current iteration 
% optimValues           values of the objective function at current iteration
% state                 current state of the optimiser

%OUTPUT
% VariableHistory       values of the optimisation variables at all iterations
% FncHistory            values of the objective function at all iterations


persistent Vhistory Fhistory
stop = false;

switch state
    case 'init'
        Vhistory = [];
        Fhistory = [];
    case 'iter'
        % Concatenate current point and objective function
        % value with history. x must be a row vector
        Fhistory = [Fhistory; optimValues.fval];
        Vhistory = [Vhistory, x(:)];
        assignin('base','VariableHistory',Vhistory);
        assignin('base','FncHistory',Fhistory);
        %if the function value changes by less than 5% in the last N iterations, stop the optimiser
%         N = 60;
%         if optimValues.iteration > N
%             Stop = 0;
%             StopCond1 = abs((Fhistory(end)-Fhistory(end-N))/Fhistory(end) *100);
%             StopCond2 = abs((Fhistory(end)-Fhistory(end-N/2))/Fhistory(end) *100);
%             if StopCond1<1 && StopCond2<1; Stop = 1; end
%             stop = Stop;
%         end
    case 'done'
        assignin('base','VariableHistory',Vhistory);
        assignin('base','FncHistory',Fhistory);
    otherwise
end
   

end


