function permutation = calcPermutation(A, B)
    [~, IA] = sort(A);
    [~, IB] = sort(B);
    permutation(IB) = IA;
end
