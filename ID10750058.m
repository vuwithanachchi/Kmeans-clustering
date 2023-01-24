
rng default;
% creates the data matrix using k-means
X = gen_kmeansdata(10750058);

disp(X)

% determine the row size & columns 
size(X,1)
size(X,2)

% getting the column amount
no_col = size(X, 2);

% mean of the column
mean_col = zeros(1, no_col);


for t = 1:no_col
    % getting the mean in column 
    mean_col(t) = mean(X(:, t));

    % printing the column mean
    fprintf('Mean of the column %d: \n', t, mean_col(t));
end
fprintf('\n');

% defining an array for standard deviation 
standard_Dev =  zeros(1, no_col);

for t = 1:no_col
    % getting the standard deviation of the each column
    standard_Dev(t) = std(X(:, t));

    % Printing the standard deviation of the each column
    fprintf('Standard deviation of the column %d: \n', t, standard_Dev(t));
end

% set up the histogram
hold on

for t = 1:no_col
    
    %  getting the histrogram of each column           
    histogram(X(:, t));  

    title(sprintf('Histogram  %d', t));    
end

hold off


covariance = cov(X);  %covariance of the matrix
correlation = corrcoef(X); %correlation of the matrix


% loop for repeating different k-means value 
K = 3:5;

% defining an array for silhouette values
value = zeros(size(K));


for t = 1:length(K)
    % getting the current value
    [c_idx, centroids] = kmeans(X, K(t));

    % getting silhouette score 
    sil = silhouette(X, c_idx);

  
    value(t) = mean(sil);

    
    figure;
    silhouette(X, c_idx);
    title(sprintf('Number of clusters = %d', K(t)));
end


fprintf('Mean silhouette value:\n');
for t = 1:length(K)
    fprintf('K = %d: %.4f\n', K(t), value(t));
end

% finding the optimal number
[~, max_index] = max(value);


fprintf('Optimal mumber of cluster %d\n', K(max_index))

% plotting the mean silhoutte score
plot(K, value, '-o')
xlabel('Number of clusters (K)')
ylabel('Mean silhouette value')

% setting the values for K
K = 3:5;


for j = K
    % Using k-means 
    [c_idx, centroids] = kmeans(X, j);

    % getting the colour map
    cmap = colormap(hsv(j));

    % ploting the clusters & the cluster centroids with colours
    figure;
    
    gscatter(X(:,1), X(:,2), c_idx, 'rbgcy','...',7)

    % ploting the clusters & the cluster centroids
    hold on;
    plot(centroids(:,1), centroids(:,2), 'kx', 'MarkerSize', 15, 'LineWidth', 3);

    title(sprintf('K_Means Clustering for the k = %d', j));

    
    leg_strings = cell(1, j);
    for t = 1:j
        leg_strings{t} = sprintf('Cluster %d', t);
    end
    legend(leg_strings);

    
    xlabel('X');% mark the x-axis 
    ylabel('Y');% mark the y-axis

    hold off;
end