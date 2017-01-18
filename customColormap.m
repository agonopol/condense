function newmap = customColormap(cmin_input, cmax_input, varargin)
    %CUSTOMCOLORMAP   Creation of a custom color map.
    %   this matlab file is designed to draw anomaly figures. the color of
    %   the colorbar is from blue to white and then to red, corresponding to 
    %   the anomaly values from negative to zero to positive, respectively. 
    %   The color white always correspondes to value zero. 
    %   
    %   You should input two values like caxis in matlab, that is the min and
    %   the max value of color values designed.  e.g. colormap(b2r(-3,5))
    %   
    %   the brightness of blue and red will change according to your setting,
    %   so that the brightness of the color corresponded to the color of his
    %   opposite number
    %   e.g. colormap(b2r(-3,6))   is from light blue to deep red
    %   e.g. colormap(b2r(-3,3))   is from deep blue to deep red
    %
    %   I'd advise you to use colorbar first to make sure the caxis' cmax and cmin.
    %   Besides, there is also another similar colorbar named 'darkb2r', in which the 
    %   color is darker.
    %
    %   by Cunjie Zhang, 2011-3-14
    %   find bugs ====> email : daisy19880411@126.com
    %   updated:  Robert Beckman help to fix the bug when start point is zero, 2015-04-08
    %   
    %   Examples:
    %   ------------------------------
    %   figure
    %   peaks;
    %   colormap(b2r(-6,8)), colorbar, title('b2r')
    %   
    
    %% check the input
    if nargin < 2 ;
       disp('input error');
       disp('input two variables, the range of caxis , for example : colormap(b2r(-3,3))');
    end
    
    % default colors: blue white red
    red_top     = [1 0 0];
    white_middle= [1 1 1];
    blue_bottom = [0 0 1];
    
    % parsing of variable argument list
    for i=1:length(varargin)-1
        if (strcmp(varargin{i}, 'presetColors'))
            switch (varargin{i+1})
                case 'ysm'
                    red_top     = [224/255 31/255 43/255];
                    white_middle= [1 1 1];
                    blue_bottom = [12/255 111/255 70/255];
                case 'ysm_gray'
                    red_top     = [224/255 31/255 43/255];
                    white_middle= [150/255 150/255 150/255];
                    blue_bottom = [12/255 111/255 70/255];
                case 'b2r'
                    red_top     = [1 0 0];
                    white_middle= [1 1 1];
                    blue_bottom = [0 0 1];
                case 'b2o_ppt'
                    red_top     = [237/255 125/255 49/255];
                    white_middle= [1 1 1];
                    blue_bottom = [91/255 155/255 213/255];
                case 'b2o_ppt_gray'
                    red_top     = [237/255 125/255 49/255];
                    white_middle= [213/255 213/255 213/255];
                    blue_bottom = [91/255 155/255 213/255];
            end
        end
        if (strcmp(varargin{i}, 'customColors'))
            red_top     = varargin{i+1};
            white_middle= varargin{i+2};
            blue_bottom = varargin{i+3};
        end
    end
    
    if cmin_input >= cmax_input
        disp('input error');
        disp('the color range must be from a smaller one to a larger one');
    end
    
    %% color interpolation 
    
    color_num = 251;   
    color_input = [blue_bottom;  white_middle;  red_top];
    oldsteps = linspace(-1, 1, size(color_input,1));
    newsteps = linspace(-1, 1, color_num);  
    
    %% Category Discussion according to the cmin and cmax input
    
    %  the color data will be remaped to color range from -max(abs(cmin_input),cmax_input)
    %  to max(abs(cmin_input),cmax_input) , and then squeeze the color data
    %  in order to make sure the blue and red color selected corresponded
    %  to their math values
    
    %  for example :
    %  if b2r(-3,6) ,the color range is from light blue to deep red , so that
    %  the light blue valued at -3 correspondes to light red valued at 3
    
    
    %% Category Discussion according to the cmin and cmax input
    % first : from negative to positive
    % then  : from positive to positive
    % last  : from negative to negative
    
    for j=1:3
       newmap_all(:,j) = min(max(transpose(interp1(oldsteps, color_input(:,j), newsteps)), 0), 1);
    end
    
    if (cmin_input < 0)  &&  (cmax_input > 0) ;  
        if abs(cmin_input) < cmax_input 
            % |--------|---------|--------------------|    
          % -cmax      cmin       0                  cmax         [cmin,cmax]
           start_point = max(round((cmin_input+cmax_input)/2/cmax_input*color_num),1);
           newmap = squeeze(newmap_all(start_point:color_num,:));
        elseif abs(cmin_input) >= cmax_input
             % |------------------|------|--------------|    
           %  cmin                0     cmax          -cmin         [cmin,cmax]   
           end_point = max(round((cmax_input-cmin_input)/2/abs(cmin_input)*color_num),1);
           newmap = squeeze(newmap_all(1:end_point,:));
        end
    elseif cmin_input >= 0
            % |-----------------|-------|-------------|    
          % -cmax               0      cmin          cmax         [cmin,cmax]
           start_point = max(round((cmin_input+cmax_input)/2/cmax_input*color_num),1);
           newmap = squeeze(newmap_all(start_point:color_num,:));
    elseif cmax_input <= 0
             % |------------|------|--------------------|    
           %  cmin         cmax    0                  -cmin         [cmin,cmax]      
           end_point = max(round((cmax_input-cmin_input)/2/abs(cmin_input)*color_num),1);
           newmap = squeeze(newmap_all(1:end_point,:));
    end
end
