% function x = load_bbfile(varargin)
% Uso:
%   x = load_bbfile(filename)
%   x = load_bbfile(filename,type)
%   x = load_bbfile(filename,type,endian)
function x = load_bbfile(varargin)
  
  switch(length(varargin))
    case 1
      filename=varargin{1};
      type='uint8';
      endian='native';
    case 2
      filename=varargin{1};
      type=varargin{2};
      endian='native';
    case 3
      filename=varargin{1};
      type=varargin{2};
      endian=varargin{3};
    otherwise
      fprintf('Error en los argumentos\n');
      fprintf('Uso : x = load_bbfile(filename)\n');
      fprintf('Uso : x = load_bbfile(filename,type)\n');
      fprintf('Uso : x = load_bbfile(filename,type,endian)\n');
      return
  end
  
  f = fopen(filename,'rb',endian);
  x = fread(f,Inf,type);
  if ~isempty(x)
    if strcmp(type,'uint8')
      x = double(x);
      x = x - 128*(x>0) + 128*(x<0);
      x = (x(1:2:end)+1j*x(2:2:end))/128.0;
    else
      x = double(x(1:2:end)+1j*x(2:2:end))/32768.0;
    end
  end 
  fclose(f);
end
