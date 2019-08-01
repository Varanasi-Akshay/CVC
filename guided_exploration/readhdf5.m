clc;clear all; close all;
addpath('/home/akshay/Desktop/CVC/HDF5/Code');
filename = 'chunked_dataset_100.hdf5';
datasetname = '/chunked';
hinfo = hdf5info(filename);
start = [200 200]; % starting point
count = [100 100]; % how many in each dimension (Chunk size)
%h5disp(filename)
data = h5read(filename,datasetname,start,count);
data_full = h5read(filename,datasetname);