function [neurons, transformed] = adjency(data)
    neurons = union(unique(data.Neuron1), unique(data.Neuron2));
    sections = unique(data.EMSection);
    transformed = zeros(length(neurons) * length(sections), length(neurons));
    for i = 1:length(sections)
        offset = (i - 1) * length(sections);
        section = data(data.EMSection == sections(i), {'Neuron1', 'Neuron2'});
        for n = 1:height(section)
            n1 = find(neurons == section(n, 'Neuron1').Neuron1);
            n2 = find(neurons == section(n, 'Neuron2').Neuron2);
            transformed(n1+offset, n2) = 1;
        end
    end
end