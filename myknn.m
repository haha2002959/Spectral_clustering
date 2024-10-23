function [ Result ] = myknn( data,k )
%计算K近邻的一个通用函数
%要求一行为一个多维变量
%返回Result为一个矩阵，第i行为到第i个数据点的欧几里得距离最短的k个数据点序号的排列，列序数越小距离越近

datasize = size(data);
dist = zeros(datasize(1),datasize(1));%距离矩阵，存放两点之间的距离
 % 计算距离矩阵
for i = 1:datasize(1)-1    
    for j = i+1:datasize(1)
        tempsum = 0;
        for ci = 1:datasize(2)
            tempsum = tempsum + (data(i,ci)-data(j,ci))*(data(i,ci)-data(j,ci));
        end
        dist(i,j) = sqrt(tempsum);
        dist(j,i) = sqrt(tempsum);
    end    
end

Result = zeros(datasize(1),k);%结果矩阵

 % 计算结果矩阵
 %datasort为距离矩阵每一行按升序排列后的结果，order为位置索引信息
 [distsort,order]=sort(dist,2,'ascend');
 Result=order(:,2:(k+1));

end
