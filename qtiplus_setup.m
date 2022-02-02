% qtiplus_setup script
function qtiplus_setup()
% call this function from within the qtiplus folder
% simply add all folders to path
addpath(genpath(pwd));
if ispc || ismac
    savepath;
    fprintf('Path set! qtiplus folder added permantently to Matlab path!\n')
else
    fprintf('Path set!\n')
end
end

