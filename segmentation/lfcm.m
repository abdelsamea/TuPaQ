function [BmuC4, center, U, ERRO] = lfcm(data, cluster_n, options)
% data=imread(filename);
% data = double(rgb2gray(data));

%The implementation of the local Fuzzy C-Means
%   [CENTER, U, OBJ_FCN] = FFCMW(DATA, N_CLUSTER) finds N_CLUSTER number of
%   clusters in the data set DATA. 

%   DATA is size M-by-N, where M is the number of
%   data points and N is the number of coordinates for each data point. The
%   coordinates for each cluster center are returned in the rows of the matrix
%   CENTER. 

%   The membership function matrix U contains the grade of membership of
%   each DATA point in each cluster. The values 0 and 1 indicate no membership
%   and full membership respectively. Grades between 0 and 1 indicate that the
%   data point has partial membership in a cluster. At each iteration, an
%   objective function is minimized to find the best location for the clusters
%   and its values are returned in OBJ_FCN.
%
%   
%       OPTIONS(1): exponent for the matrix U             (default: 2.0)
%       OPTIONS(2): maximum number of iterations          (default: 100)
%       OPTIONS(3): minimum amount of improvement         (default: 1e-5)
%       OPTIONS(4): info display during iteration         (default: 1)
%   The clustering process stops when the maximum number of iterations
%   is reached, or when the objective function improvement between two
%   consecutive iterations is less than the minimum amount of improvement
%   specified. 

if nargin ~= 2 && nargin ~= 3,
	error('Too many or too few input arguments!');
end
%% local data 
%data= imgaussfilt(data,1);
sigma=1.;
KernelSize = 2*round(2*sigma)+1;
K = fspecial('gaussian', KernelSize, sigma); 
data = conv2(data,K,'same');

[a1 a2]=size(data);
data = data(:);
data_n = size(data, 1);
%%

% Change the following to set default options
default_options = [2;	% exponent for the partition matrix U
		200;	% max. number of iteration
		1e-5;	% min. amount of improvement
		1];	% info display during iteration 

if nargin == 2,
	options = default_options;
else
	% If "options" is not fully specified, pad it with default values.
	if length(options) < 4,
		tmp = default_options;
		tmp(1:length(options)) = options;
		options = tmp;
	end
	% If some entries of "options" are nan's, replace them with defaults.
	nan_index = find(isnan(options)==1);
	options(nan_index) = default_options(nan_index);
	if options(1) <= 1,
		error('The exponent should be greater than 1!');
	end
end

expo = options(1);		    % Exponent for U
max_iter = options(2);		% Max. iteration
min_impro = options(3);		% Min. improvement
display = options(4);		% Display info or not

obj_fcn = zeros(max_iter, 1);	% Array for objective function

U = rand(cluster_n, data_n);
col_sum = sum(U);
U = U./col_sum(ones(cluster_n, 1), :);
sumDdotD2 = sum(data.*data,2)';

% Main loop
ERRO =0;
for i = 1:max_iter,
    mf = U.^expo;       % MF matrix after exponential modification    
    center = bsxfun(@rdivide,mf*data,sum(mf,2));
    dist = bsxfun(@plus,sumDdotD2,sum(center.*center,2))-2*(center*data');        
    aux = dist.*mf;
    obj_fcn(i) = sum(aux(:));  % objective function
    tmp = dist.^(-1/(expo-1));
    U = bsxfun(@rdivide, tmp, sum(tmp,1));

	%if display, 
	%	fprintf('Iteration count = %d, obj. fcn = %f\n', i, obj_fcn(i));
	%end
	% check terminating condition
	
    
   if i > 1,
       ERRO = abs(obj_fcn(i) - obj_fcn(i-1));
		if ERRO < min_impro, break; end,
   end
end

iter_n = i;	% Actual number of iterations 
obj_fcn(iter_n+1:max_iter) = [];

 %[q , indexx]=min(sum(U,2));
     [q , indexx]=min(center);
     BmuC4=reshape(U(indexx,:),a1,a2,1);
     %BmuC4(isnan(BmuC4)) = 0;
% BmuC1=reshape(U(1,:),a1,a2,1);
% BmuC2=reshape(U(2,:),a1,a2,1);
% BmuC3=reshape(U(3,:),a1,a2,1);
% figure(500), imshow(BmuC1,[]);
% figure(501), imshow(BmuC2,[]);
% figure(502), imshow(BmuC3,[]);
% figure(503), imshow(BmuC4,[]);
%  if sum(U(1,:)) > sum(U(2,:))
%      BmuC1=reshape(U(1,:),a1,a2,1);
%  else
%      BmuC1=reshape(U(2,:),a1,a2,1);
%  end
% figure(1001), imshow(BmuC1,[]);

end