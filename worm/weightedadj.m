function [neurons, matrix] = weightedadj(data)
    neurons = union(unique(data.Neuron1), unique(data.Neuron2));
    sections = unique(data.EMSection);
    matrix = zeros(length(neurons), length(neurons));
    for i = 1:length(sections)
        section = data(data.EMSection == sections(i), {'Neuron1', 'Neuron2'});
        for n = 1:height(section)
            n1 = find(neurons == section(n, 'Neuron1').Neuron1);
            n2 = find(neurons == section(n, 'Neuron2').Neuron2);
            matrix(n1,n2) = matrix(n1,n2) + 1;
            matrix(n2,n1) = matrix(n2,n1) + 1;
        end
    end
end