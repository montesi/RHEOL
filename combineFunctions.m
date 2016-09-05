function fh = combineFunctions(funcHandles)
    %# return a function handle
    fh = @myGeneralFunction;

    %# nested function. Implements the formula:
    %# f(x) = cos( f1(x) + f2(x) + ... + fN(x) )
    %# where f1,..,fN are the passed function handles 
    function y = myGeneralFunction(x)
        %# evaluate all functions on the input x
        y = cellfun(@(fcn) fcn(x), funcHandles);

        %# apply cos(.) to the sum of all the previous results
        %# (you can make this any formula you want)
        y =  sum(y) ;
    end
end