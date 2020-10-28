%% Random parameter generator for data validation
% Generate a text file with 36 random number from 0 to 1

path='C:\Data\';
fileName='fileName';

n=36;
M=rand([n 1]);
dlmwrite([path,fileName,'.txt'],M,'delimiter','\t')