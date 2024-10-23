function [index_cluster,cluster] = kmeans_func(data,cluster_num)
[m,n]=size(data);
cluster=data(randperm(m,cluster_num),:);%��m���������ѡ��cluster_num������Ϊ��ʼ�������ĵ�
epoch_max=1000;%����������
therad_lim=0.001;%���ı仯��ֵ,����������ε����ľ������ĵ����֮��С�ڸ���ֵ��������ࡣ
epoch_num=0;
while(epoch_num<epoch_max)
    epoch_num=epoch_num+1;
    % distance1�洢ÿ���㵽���������ĵ�ŷ�Ͼ���
    for i=1:cluster_num
        distance=(data-repmat(cluster(i,:),m,1)).^2;
        distance1(:,i)=sqrt(sum(distance'));
    end
    %��ʣ���ÿ���㶼�������������ľ����������ڵ����
    %index_cluster���ڴ洢ÿһ�е���Сֵ���ڵ��к���������Ϊÿ�������ر�ǩ
    [~,index_cluster]=min(distance1');%index_clusterȡֵ��Χ1~cluster_num
    % cluster_new�洢�µľ�������
    for j=1:cluster_num
        cluster_new(j,:)=mean(data(find(index_cluster==j),:));
    end
    %����µľ������ĺ���һ�ֵľ������ľ���ʹ������ı仯��ֵtherad_lim�����¾������ģ������㷨����
    if (sqrt(sum((cluster_new-cluster).^2))>therad_lim)
        cluster=cluster_new;
    else
        break;
    end
end
end