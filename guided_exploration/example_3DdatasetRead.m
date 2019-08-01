function data_out = example_3DdatasetRead
% Example to read part of a 3D dataset

filename = 'threeD.h5';
fileID = H5F.open(filename, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');

datasetName = '/dataset1';
datasetID = H5D.open(fileID, datasetName);

% Get dataspace
dataspaceID = H5D.get_space(datasetID);
[rank dims] = H5S.get_simple_extent_dims(dataspaceID);

% Select hyperslab of data 
dimsRequired = [dims(1) dims(2) dims(3)/5]; % L x M x N/5
start = [0 0 0];
stride = [];
count = dimsRequired;
block = [];
H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', ...
                     start, stride, count, block);
                 
% Define the memory dataspace.
memspaceID = H5S.create_simple(rank, count, []);

% Read the subsetted data
data_out = H5D.read(datasetID, 'H5ML_DEFAULT', ...
                    memspaceID, dataspaceID, 'H5P_DEFAULT');
                

