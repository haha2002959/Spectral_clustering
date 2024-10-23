function [ Result ] = myknn( data,k )
%����K���ڵ�һ��ͨ�ú���
%Ҫ��һ��Ϊһ����ά����
%����ResultΪһ�����󣬵�i��Ϊ����i�����ݵ��ŷ����þ�����̵�k�����ݵ���ŵ����У�������ԽС����Խ��

datasize = size(data);
dist = zeros(datasize(1),datasize(1));%������󣬴������֮��ľ���
 % ����������
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

Result = zeros(datasize(1),k);%�������

 % ����������
 %datasortΪ�������ÿһ�а��������к�Ľ����orderΪλ��������Ϣ
 [distsort,order]=sort(dist,2,'ascend');
 Result=order(:,2:(k+1));

end
