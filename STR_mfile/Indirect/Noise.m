clc
clear all
close all
t=500;
varN=0.01;% noise variance
r=random('normal',0,varN,t,1);
noise_minor=[zeros(1,3) r']';