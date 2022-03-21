function [kdotuR,kdotuI] = dotProduct(kx,ky,kz,uhatR,uhatI,vhatR,vhatI,whatR,whatI)
    kdotuR = kx.*uhatR + ky.*vhatR + kz.*whatR;
    kdotuI = kx.*uhatI + ky.*vhatI + kz.*whatI;