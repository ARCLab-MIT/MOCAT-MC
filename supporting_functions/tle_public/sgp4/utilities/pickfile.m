function filename = pickfile ( filter , title, indir, permission )
% The purpose of this function simplifies the user interface to uigetfile.
% If filter is undefined, '*.*' will be used as the default.
% The filename returned is [pname fname] returned by uigetfile.
% 

persistent dirname

if ~exist('filter','var')
	filter = '*.*';
end
if ~exist('title','var')
	title = 'Select File to Open';
end
if ~exist('indir','var')
   use_indir = 0;
else
   if isempty(indir)
      use_indir = 0;
   else
      if (exist(indir)==7) use_indir = 1; else use_indir = 0; end
   end   
end
if ~exist('permission','var')
   permission = 'r';
end

current_dir = pwd;

if use_indir
   cd(indir);
else
   if ~isempty(dirname)
      if (exist(dirname) == 7) cd(dirname); end;
   end
end

if strcmpi(permission,'w')
   [fname, pname] = uiputfile(filter,title);
else   
   [fname, pname] = uigetfile(filter,title);
end
filename = [pname fname];

if (pname == 0)
   filename = 0;
   if ~use_indir dirname = current_dir; end;
else
   if ~use_indir dirname = pname; end;
end

cd(current_dir);
   
return

