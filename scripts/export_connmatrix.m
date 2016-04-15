function export_connmatrix(name,m,msize,popsize,offset,batch_index,batch_size,output_dir)

  range = offset+(batch_index:batch_size:popsize);
  chunkSize=10000000;
  chunkStart=1; chunkEnd = chunkSize+chunkStart;
  filename = sprintf('%s/%s.%d.dat', output_dir, name, batch_index);
  fid = fopen (filename,'w');
  fclose(fid);
  while (chunkStart < msize)
    if chunkEnd < msize
       data = double(m(ismember(m(chunkStart:chunkEnd,2),range),:));
    else
       data = double(m(ismember(m(chunkStart:msize,2),range),:));
    end
    fprintf('export_connmatrix: name = %s offset = %d popsize = %d size(data) = %d batch_index = %d\n', name, offset, popsize, size(data,1), batch_index)
    if (not (isempty(data)))
      save('-ascii', '-append', filename, 'data');
    end
    chunkStart = chunkEnd+1;
    chunkEnd = chunkSize+chunkStart;
  end
end
