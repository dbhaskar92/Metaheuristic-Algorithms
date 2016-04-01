% Parse ORLIB BQP test problems

filename = 'bqp500.txt';

fileID = fopen(filename);

totalProbs = textscan(fileID, '%d', 1);

for i = 1 : totalProbs{1,1}
   
    Qinfo = textscan(fileID, '%d %d', 1);
    n = Qinfo{1,1};
    rows2parse = Qinfo{1,2};
    
    Q = zeros(n,n);
    
    Matrix = textscan(fileID, '%d %d %d', rows2parse);
    indexi = Matrix{1,1};
    indexj = Matrix{1,2};
    values = Matrix{1,3};
    
    for j = 1 : rows2parse
        
        Q(indexi(j), indexj(j)) = values(j);
        Q(indexj(j), indexi(j)) = values(j);
     
    end
    
    % write the matrix to file
    writefilename = strcat('m', int2str(i), filename);
    dlmwrite(writefilename, Q);
        
end
