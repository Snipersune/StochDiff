function newDW = coarsen(dW,k)
    % Coarsen the discretized dW with a factor k
    numProc = size(dW,1);
    N = size(dW,2);

    newN = cast(floor(N/k), "int64");
    newDW = zeros(numProc,newN);
    for i = 1:newN
        indexVec = (k*(i-1)+1) : (k*(i-1)+k);
        newDW(:,i) = sum(dW(:,indexVec),2);
    end
end