function y = parameterized_objective(s,x)

% This function calls the accounting model, cashflows,
% and is called by the ObjectiveFunction fuction
% handle, which is itself called within ga.

[cashflows_outputs] = cashflows(s,x);

y = cashflows_outputs.y;

end

% I am not sure why this function is necessary, instead
% of calling cashflows directly, like this:
% ObjectiveFunction = @(x) cashflows(s,x)
% But it is a trivial change that would bring close to
% no benefit.