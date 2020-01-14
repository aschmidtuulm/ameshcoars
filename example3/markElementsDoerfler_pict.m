function marked = markElementsDoerfler_pict(eta,theta)

[S,I] = sort(eta,'descend'); %sortiert absteigend und speichert Indices in I
i = find(cumsum(S) > theta*sum(eta),1);
marked = I(1:i);

end