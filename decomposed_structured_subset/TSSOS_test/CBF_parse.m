function model = CBF_parse(filepath)
  %https://github.com/HFriberg/micoso-solver/blob/master/%40CBFdata/parse.m
  %slight modifications
  model = struct;
  model.ver = 0;

  f = fopen(filepath, 'rt');
  if (f == -1)
    error(['fopen failed: ', filepath])
  end
  
  TXT = fscanf(f, '%c');
  NEWLINES = find(TXT == 10);
  ST = [1, NEWLINES(1:end-1)+1];
  EN = NEWLINES(1:end)-1;
  
  DATALINES = (ST<=EN & TXT(ST)~='#');
  ST = ST(DATALINES);
  EN = EN(DATALINES);
  
  fclose(f);

  li = 0;
  while li < length(ST)
    [li,line] = next(TXT,ST,EN,li,1);
    
    if (model.ver == 0)  
      if strcmp(line, 'VER')
        [li,line] = next(TXT,ST,EN,li,1);
        model.ver = sscanf(line, '%i');
        %if (model.ver ~= 1)
        %  error('The version of the file format is not supported.')
        %end
        %attempting to do for version 3
      else
        error('First keyword should be VER.')
      end
      
    else
      
      switch(line)
        case 'OBJSENSE'
          [li,line] = next(TXT,ST,EN,li,1);
          model.objsense = sscanf(line, '%s');
          
        case 'VAR'
          [li,line] = next(TXT,ST,EN,li,1);
          buf = sscanf(line, '%u %u');
          model.varnum = buf(1); 
          varstacknum = buf(2);
          
          [li,line] = next(TXT,ST,EN,li,varstacknum);
          buf = textscan(line, '%s %u');
          model.varstackdomain = buf{1};
          model.varstackdim = double(buf{2});
          
        case 'INT'
          [li,line] = next(TXT,ST,EN,li,1);
          intvarnum = sscanf(line, '%u');
          model.intvar = zeros(intvarnum, 1);
          
          [li,line] = next(TXT,ST,EN,li,intvarnum);
          buf = textscan(line, '%u');
          model.intvar = double(buf{1}) + 1;
          
        case 'CON'
          [li,line] = next(TXT,ST,EN,li,1);
          buf = sscanf(line, '%u %u');
          model.mapnum = buf(1); 
          mapstacknum = buf(2);
          
          [li,line] = next(TXT,ST,EN,li,mapstacknum);
          buf = textscan(line, '%s %u');
          model.mapstackdomain = buf{1};
          model.mapstackdim = double(buf{2});
          
        case 'PSDVAR'
          error('NO SUPPORT')
          
        case 'PSDCON'
          error('NO SUPPORT')

        case 'OBJFCOORD'
          error('NO SUPPORT')
          
        case 'OBJACOORD'
          [li,line] = next(TXT,ST,EN,li,1);
          objannz = sscanf(line, '%u');
          
          [li,line] = next(TXT,ST,EN,li,objannz);
          buf = textscan(line, '%u %f');
          
          model.c = sparse(double(buf{1})+1, 1, buf{2}, ...
                          model.varnum, 1, objannz);
          
        case 'OBJBCOORD'
          [li,line] = next(TXT,ST,EN,li,1);
          model.c0 = sscanf(line, '%f');
          
        case 'FCOORD'
          error('NO SUPPORT')
          
        case 'ACOORD'
          [li,line] = next(TXT,ST,EN,li,1);
          annz = sscanf(line, '%u');
          
          [li,line] = next(TXT,ST,EN,li,annz);
          buf = textscan(line, '%u %u %f');
          
          model.A = sparse(double(buf{1})+1, double(buf{2})+1, buf{3}, ... 
                          model.mapnum, model.varnum, annz);
          
        case 'BCOORD'
          [li,line] = next(TXT,ST,EN,li,1);
          bnnz = sscanf(line, '%u');
          
          [li,line] = next(TXT,ST,EN,li,bnnz);
          buf = textscan(line, '%u %f');
          
          model.b = sparse(double(buf{1})+1, 1, buf{2}, ...
                          model.mapnum, 1, bnnz);
          
        case 'HCOORD'
          error('NO SUPPORT')
          
        case 'DCOORD'
          error('NO SUPPORT')
          
        otherwise
          error(['Keyword "', line, '" not recognized!'])
      end
    end
  end
  
  if (isempty(model.c) && model.varnum >= 1)
    model.c = sparse(model.varnum, 1);
  end
  
  if strcmp(model.objsense, "MAX")
      model.c = -model.c;
  end
  
  if (isempty(model.b) && model.mapnum >= 1)
    model.b = sparse(model.mapnum, 1);
  end
  
  
end

function [li,line] = next(TXT,ST,EN,li,num)
  line = TXT(ST(li+1):EN(li+num));
  li = li + num;
end
