function [MINU] = LoadNeuFeature(fname)

data = textread(fname);
% sp_num = data(3,1)+data(4,1);
% MINU = data((6+sp_num):end,1:3);
MINU = data(1:end,1:3);
MINU(:,3) = - MINU(:,3);

end

