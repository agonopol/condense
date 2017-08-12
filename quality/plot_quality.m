function plot_quality(data, fn, varargin) 
   output = 'quality.png';
   hold on;
   for i=1:2:length(varargin)-1
        if (strcmp(varargin{i}, 'output'))
            output = varargin{i+1};
        else 
            label = varargin{i};
            assigments = varargin{i+1};
            x = arrayfun(@(i) max(assigments(i, :)), 1:size(assigments,1));
            y = qscore(data, assigments, 'quality', fn);
            plot (x, y, 'DisplayName' , label);
        end
   end
   hold off;
   legend ( 'show' );
   
   saveas(gcf, output);
   close all;
end