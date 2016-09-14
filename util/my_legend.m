function my_legend(cellStrings)
%MY_LEGEND Summary of this function goes here
%   Detailed explanation goes here

h_legend=legend(cellStrings);
set(h_legend,'FontWeight', 'bold' , 'FontSize',20, 'FontName', 'Arial');


end

