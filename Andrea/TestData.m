function TestData()

   
    %data1 = readmatrix("HIVIn(Matlab).csv");
    data2 = readmatrix("LeukB4(Matlab).csv")


    %setup(data1, "step", "HIV")

    data2 = data2(sum(isnan(data2),2)==0,:);
    setup(data2, "step", "LeukB4")

    %stepminer test
    %index is mean 1, mean 2, breakpoint 1, breakpoint 2
   %  h = height(bigdata)
   %  n = width(bigdata)
   % 
   % kmeansData = zeros(h, n)
   % BASCdata = zeros(h, n)
   % stepminerData = zeros(h, n)
   % %BASCBdata = zeros(h, n)
   % 
   %  for i = 1:h
   % 
   %      %read array 
   %      array = bigdata(i,:)
   % 
   %      %alpha 0.05 for test
   %      kmeansData(i,:) = kmeans(array,2)
   %      BASCdata(i,:) = BASCA(array)
   %      [type(i), stepminerData(i,:)] = stepminer(array,0.05)
   %      %BASCBdata(i,:) = BASCB(array)
   %  end
   % 
   % 
   %  writematrix(kmeansData, "kmeans.txt")
   %  writematrix(BASCdata, "BASCA.txt")
   %  %writematrix(BASCBdata, "BASCB.txt")
   %  writematrix(stepminerData, "stepminer.txt")
   %  writematrix(transpose(type), "steptype.txt")
   % 


end


function setup(x, fun, name)

    h = height(x);
    n = width(x);

    
    for i = 1:h
   
         %read array 
         array = x(i,:)
        
         splineGraph(array,name,fun,i,n)
     end

end



function splineGraph(array,name,fun,index,n)
    
   %spline

    pp = spline(1:n, array);
    yy1 = ppval(pp,1:n/10:n);
    yy2 = ppval(pp,1:n/20:n);
    yy4 = ppval(pp,1:n/40:n);
    yy8 = ppval(pp,1:n/80:n);

    if fun == "kmeans"

        [~,t1] = kmeans(yy1,2);
        [~,t2] = kmeans(yy2,2);
        [~,t4] = kmeans(yy4,2);
        [~,t8] = kmeans(yy8,2);

    elseif fun == "BASC"

        [~,t1] = BASCA(yy1);
        [~,t2] = BASCA(yy2);
        [~,t4] = BASCA(yy4);
        [~,t8] = BASCA(yy8);

    elseif fun == "step"

        [~,t1] = stepminer(yy1,0.05);
        [~,t2] = stepminer(yy2,0.05);
        [~,t4] = stepminer(yy4,0.05);
        [~,t8] = stepminer(yy8,0.05);
    end

    plot(yy8)
    y1 = yline(t1,'-r', '10');
    y1.LabelHorizontalAlignment ='center';
    y1.LabelVerticalAlignment='top';
    y2 = yline(t2,'--y','20');
    y2.LabelHorizontalAlignment='left';
    y4 = yline(t4,':g','40');
    y4.LabelHorizontalAlignment = 'right';
    y8 = yline(t8,'-.b','80');
    y8.LabelHorizontalAlignment ='center';
    y8.LabelVerticalAlignment='bottom';
    
    f = gca;
    filename = strcat('img/',name,'/', fun,'/',num2str(index,'%02d'),".png")
    exportgraphics(f,filename)
end