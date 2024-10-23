clc, clear
X = load('data_moon.csv');
plot(X(:,1), X(:,2), 'k.');
%% 数据初始化
[n,~] = size(X); % 获取样本数n
W0 = zeros(n,n); % 全连接邻接矩阵
W = zeros(n,n);  % 经过k近邻提取后的亲和度(边权重)矩阵
D = zeros(n,n);  % 度矩阵,亲和矩阵 每行的和 被赋值到主对角线
Knear = 10;      % 前k个近邻
k_eigvec = 2;    % 前k小特征值对应特征向量
sigma = 0.8;     % 高斯核中的标准差
classCnt = 2;    % 聚类的个数
%% 1.图构造
% 1.1 构造一个全连接的邻接矩阵W0,使用高斯核作为距离度量,点对距离越近值越大
for i = 1:n
    for j = 1:n
        if i ~= j
            W0(i,j) = exp( -(norm(X(i,:) - X(j,:))^2)/(2*sigma.^2) ); % 第i个样本和第j个样本的距离
        end
    end
end
% 1.2 构造亲和度矩阵W，使用k近邻对以上全连接邻接矩阵W0的每行进行降序排序，提取前k个值作为边(i,j)的权值(亲和度)
for i = 1:n
    [~,maxKIndex] = sort(W0(i,:),'descend');  % 将第i个点到其余n-1个点的高斯距离降序排序，取出降序后的索引maxKIndex
    W(i,maxKIndex(1 : Knear)) = W0(i,maxKIndex(1 : Knear)); % 将前Knear大个样本作为第i个样本的邻居，其余点与该点则无边相连
end
% 1.3 使不对称的W矩阵变为实对称矩阵
W = (W' + W)/2;
%%  2.计算度矩阵D
for i = 1:n
    D(i,i) = sum(W(i,:));
end
%%  3.构造拉普拉斯矩阵，并进行对称型归一化
L = D - W;
Lsym = D^(-0.5) * L * D^(-0.5);
%%  4.计算Lsym的前k_eigvec小的特征值和特征向量
[eigVecCol, eigValueDig] = eig(Lsym);  % 获取Lsym的所有特征向量(列向量返回)和特征值(对角矩阵返回)
eigValue = eigValueDig * ones(n,1); % 将特征值对角矩阵转换为一个列向量，方便排序
[~,minKIndex] = sort(eigValue, 'ascend');    % 特征值升序排序
%%  5.将前k小特征向量组合成U矩阵，并按行归一化得到T矩阵
U = eigVecCol(:,minKIndex(2 : k_eigvec + 1));   % 取出最前k_eigvec小个特征值(不包含0特征值)对应的索引,组合对应特征向量构建矩阵U
T = zeros(n,k_eigvec);
for i = 1:n
    for j = 1:k_eigvec
        T(i,j) = U(i,j)/norm(U(i));
    end
end
%%  6.使用k-means做最后的聚类
[label, cluster] = kmeans_func(T,classCnt);
%%  7.输出
%画出聚类结果
% figure(1);
% cmap = colormap;
% for i = 1:classCnt
%     tmp_data = X(label==i,:);
%     ic = int8((i*64)/(classCnt*1));
%     col = cmap(ic,:);
%     plot(tmp_data(:,1),tmp_data(:,2),'o','MarkerSize',2,'MarkerFaceColor',col,'MarkerEdgeColor',col);
%  %   set(gca,'XTick',[0:5:15]) %改变x轴坐标间隔显示 这里间隔为2
%     hold on;
% end

x_1 = X(label == 1,:);
x_2 = X(label == 2,:);
%x_3 = X(label == 3,:);
plot(x_1(:,1), x_1(:,2), 'r.'); hold on; plot(x_2(:,1), x_2(:,2), 'b.');
%plot(x_3(:,1), x_3(:,2), 'y.')
title('谱聚类效果,Knear=10,sigma=0.8')
