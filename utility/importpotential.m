% parse the three body txt file for the relevant threebody parameters

function [mab,mbc,D] = importpotential(vinput)

parameters = importdata(vinput);

mab = parameters(1:4,1);
mbc = parameters(5:8,1);
D = parameters(9:end,:);    