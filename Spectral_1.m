clc, clear
X = load('data_moon.csv');
plot(X(:,1), X(:,2), 'k.');
%% ���ݳ�ʼ��
[n,~] = size(X); % ��ȡ������n
W0 = zeros(n,n); % ȫ�����ڽӾ���
W = zeros(n,n);  % ����k������ȡ����׺Ͷ�(��Ȩ��)����
D = zeros(n,n);  % �Ⱦ���,�׺;��� ÿ�еĺ� ����ֵ�����Խ���
Knear = 10;      % ǰk������
k_eigvec = 2;    % ǰkС����ֵ��Ӧ��������
sigma = 0.8;     % ��˹���еı�׼��
classCnt = 2;    % ����ĸ���
%% 1.ͼ����
% 1.1 ����һ��ȫ���ӵ��ڽӾ���W0,ʹ�ø�˹����Ϊ�������,��Ծ���Խ��ֵԽ��
for i = 1:n
    for j = 1:n
        if i ~= j
            W0(i,j) = exp( -(norm(X(i,:) - X(j,:))^2)/(2*sigma.^2) ); % ��i�������͵�j�������ľ���
        end
    end
end
% 1.2 �����׺ͶȾ���W��ʹ��k���ڶ�����ȫ�����ڽӾ���W0��ÿ�н��н���������ȡǰk��ֵ��Ϊ��(i,j)��Ȩֵ(�׺Ͷ�)
for i = 1:n
    [~,maxKIndex] = sort(W0(i,:),'descend');  % ����i���㵽����n-1����ĸ�˹���뽵������ȡ������������maxKIndex
    W(i,maxKIndex(1 : Knear)) = W0(i,maxKIndex(1 : Knear)); % ��ǰKnear���������Ϊ��i���������ھӣ��������õ����ޱ�����
end
% 1.3 ʹ���ԳƵ�W�����Ϊʵ�Գƾ���
W = (W' + W)/2;
%%  2.����Ⱦ���D
for i = 1:n
    D(i,i) = sum(W(i,:));
end
%%  3.����������˹���󣬲����жԳ��͹�һ��
L = D - W;
Lsym = D^(-0.5) * L * D^(-0.5);
%%  4.����Lsym��ǰk_eigvecС������ֵ����������
[eigVecCol, eigValueDig] = eig(Lsym);  % ��ȡLsym��������������(����������)������ֵ(�ԽǾ��󷵻�)
eigValue = eigValueDig * ones(n,1); % ������ֵ�ԽǾ���ת��Ϊһ������������������
[~,minKIndex] = sort(eigValue, 'ascend');    % ����ֵ��������
%%  5.��ǰkС����������ϳ�U���󣬲����й�һ���õ�T����
U = eigVecCol(:,minKIndex(2 : k_eigvec + 1));   % ȡ����ǰk_eigvecС������ֵ(������0����ֵ)��Ӧ������,��϶�Ӧ����������������U
T = zeros(n,k_eigvec);
for i = 1:n
    for j = 1:k_eigvec
        T(i,j) = U(i,j)/norm(U(i));
    end
end
%%  6.ʹ��k-means�����ľ���
[label, cluster] = kmeans_func(T,classCnt);
%%  7.���
%����������
% figure(1);
% cmap = colormap;
% for i = 1:classCnt
%     tmp_data = X(label==i,:);
%     ic = int8((i*64)/(classCnt*1));
%     col = cmap(ic,:);
%     plot(tmp_data(:,1),tmp_data(:,2),'o','MarkerSize',2,'MarkerFaceColor',col,'MarkerEdgeColor',col);
%  %   set(gca,'XTick',[0:5:15]) %�ı�x����������ʾ ������Ϊ2
%     hold on;
% end

x_1 = X(label == 1,:);
x_2 = X(label == 2,:);
%x_3 = X(label == 3,:);
plot(x_1(:,1), x_1(:,2), 'r.'); hold on; plot(x_2(:,1), x_2(:,2), 'b.');
%plot(x_3(:,1), x_3(:,2), 'y.')
title('�׾���Ч��,Knear=10,sigma=0.8')
