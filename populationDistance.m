function distance = populationDistance(populations, method)
   if size(populations, 1) == 2
        distance = method(populations{1}, populations{2});
   else
       distance = arrayfun(@(pair) method(pair{i1}, pair{2}), combnk(populations, 2));
       distance = mean(distance);
   end
end