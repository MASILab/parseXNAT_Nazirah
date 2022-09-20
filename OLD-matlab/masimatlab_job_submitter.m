function [status,result] = masimatlab_job_submitter(cmd,varargin)
%Masi Matlab Job Submitter
%For use with the masi cluster
%
%Wrapper for python jobsubmitter to simplify submitting cluster jobs of
%matlab functions
%
%USAGE: result = masimatlab_job_submitter(cmd, [opts], [arg1], [arg2], ...)
%    INPUTS:
%     -cmd: string/fcn handle containing the command to submit
%          NOTE: If you are passing in strings, escape your quotes ie,
%          ['command(''',string,''')']
%
%     -opts: a structure of options. Possible options are:
%       -window: bool, whether or not to use xvfb-run (Does your command
%                require a window server?), default=false
%       -working_dir: string, directory to start command, default pwd
%       -log_dir: string, directory for log files, default pwd/logs
%       -combine_output: bool, combine stdout and stderr into one
%       file,default=false
%       -name: string, name for job, default=current time
%       -bashrc: bool, execute bashrc? default=false
%       -dont_submit: bool, default=false
%       -queue: string, queue to submit to, default clusterjob
%       -pbs_file: string, file to save pbs, default=don't save
%       -print_host_info: bool,print out machine info, default=false
%       -print_time: bool, print out duration at completion of job
%       -memory: string,memory required for job, default=4G
%       -job_hold: name or jid of jobs to wait for (string)
%       -email: email address to receive notifications
%       -args: cell array of arguments, if none are given with a function
%       handle it is assumed to be a script, can also be passed in as
%       arguments after opts.
%
%    OUTPUTS:
%       -status: status returned by job_submitter of whether completion was
%       successful
%       -return: the text returned by job_submitter, contains job number
%
%Rob Harrigan
%October 31, 2014 - Original creation
%November 19,2015 - updated to work with function handles and input
%arguments
%TODO: type check input options (ie ischar(email_addr)

%%Check inputs...lots of inputs
narginchk(1,inf);
if nargin>1
    opts = varargin{1};
end

if ~isa(opts,'struct')
    error('Options must be passed in a structure');
end

if isa(cmd, 'function_handle')
    %Check for arguments
    if nargin>2
        %arguments have been passed in not as a structure
        args = varargin(2:end);
    end
    if isfield(opts,'args')
        %arguments passed in as a structure
        if exist(args,'var')~=0
            error('Arguments must be passed in as either a structure or from the command');
        end
        if iscell(args)
            args = opts.args;
        else
            args = {opts.args};
        end
    end
    %Now build the command based on the arguments
    if exist('args','var')==0
        %It is a script, just run the command
        cmd = func2str(cmd);
    else
        cmd = build_cmd_string(cmd,args);
    end
    
end

if isfield(opts,'window')
    window = opts.window;
    opts = rmfield(opts,'window');
else
    window = false;
end

if isfield(opts,'working_dir')
    working_dir = opts.working_dir;
    opts = rmfield(opts,'working_dir');
else
    working_dir = pwd;
end
if ~isdir(working_dir)
    warning('WARNING: working directory does not exist, creating');
    mkdir(working_dir)
end

if isfield(opts,'log_dir')
    log_dir = opts.log_dir;
    opts = rmfield(opts,'log_dir');
else
    log_dir = fullfile(pwd,'logs');
end
if ~isdir(log_dir)
    mkdir(log_dir)
end

if isfield(opts,'combine_output')
    combine_output = opts.combine_output;
    opts = rmfield(opts,'combine_output');
else
    combine_output = false;
end

if isfield(opts,'name')
    name = opts.name;
    opts = rmfield(opts,'name');
else
    name = '';
end

if isfield(opts,'bashrc')
    bashrc = opts.bashrc;
    opts = rmfield(opts,'bashrc');
else
    bashrc = false;
end

if isfield(opts,'dont_submit')
    dont_submit = opts.dont_submit;
    opts = rmfield(opts,'dont_submit');
else
    dont_submit=false;
end

if isfield(opts,'queue')
    queue = opts.queue;
    opts = rmfield(opts,'queue');
else
    queue = 'clusterjob';
end

if isfield(opts,'pbs_file')
    pbs_file = opts.pbs_file;
    opts = rmfield(opts,'pbs_file');
else
    pbs_file = '';
end

if isfield(opts,'print_host_info')
    print_host_info = opts.print_host_info;
    opts = rmfield(opts,'print_host_info');
else
    print_host_info=false;
end

if isfield(opts,'print_time')
    print_time = opts.print_time;
    opts = rmfield(opts,'print_time');
else
    print_time = false;
end

if isfield(opts,'memory')
    memory = opts.memory;
    opts = rmfield(opts,'memory');
else
    memory='4G';
end

if isfield(opts,'email')
    sendEmail = true;
    email_addr = opts.email;
    opts = rmfield(opts,'email');
else
    sendEmail=false;
end

if isfield(opts,'job_hold')
    hold_job=true;
    hold_jid=opts.job_hold;
    opts = rmfield(opts,'job_hold');
else
    hold_job=false;
end

%% Build Command
%       -working_dir: string, directory to start command, default pwd
%       -log_dir: string, directory for log files, default pwd/logs
%       -combine_output: bool, combine stdout and stderr into one
%       file,default=false
%       -name: string, name for job, default=current time
%       -bashrc: bool, execute bashrc? default=false
%       -dont_submit: bool, default=false
%       -queue: string, queue to submit to, default clusterjob
%       -pbs_file: string, file to save pbs, default=don't save
%       -print_host_info: bool,print out machine info, default=false
%       -memory: string,memory required, default=4G
job_submitter_cmd = '/share/bin/job_submitter';

command_cell = cell(0);
%Call job submitter
command_cell{    1} = job_submitter_cmd;
%Command
command_cell{end+1} = ' --command ';
command_cell{end+1} = '"';
if window
    command_cell{end+1} = 'xvfb-run -a --server-args=\"-screen 0 1600x1280x24 -ac -extension GLX\"';
end
command_cell{end+1} = ' matlab -nodesktop -nosplash -singleCompThread';
command_cell{end+1} = ['-r \"',cmd,'\""'];
%Starting dir
command_cell{end+1} = ' --starting-dir ';
command_cell{end+1} = working_dir;
%log dir
command_cell{end+1} = ' --log-dir ';
command_cell{end+1} = log_dir;
%combine_output
if combine_output
    command_cell{end+1} = ' --combine-output ';
end
if ~isempty(name)
    command_cell{end+1} = ' --name ';
    command_cell{end+1} = name;
end
if bashrc
    command_cell{end+1} = ' --execute-bashrc ';
end
if dont_submit
    command_cell{end+1} = ' --dont-submit ';
end
command_cell{end+1} = ' --queue ';
command_cell{end+1} = queue;
if ~isempty(pbs_file)
    command_cell{end+1} = ' --pbs-file ';
    command_cell{end+1} = pbs_file;
end
if print_host_info
    command_cell{end+1} = ' --print-host-info ';
end
if print_time
    command_cell{end+1} = ' --print-time ';
end

if hold_job
    command_cell{end+1} = ' --hold-job ';
    command_cell{end+1} = hold_jid;
end

if sendEmail
    command_cell{end+1} = ' --email ';
    command_cell{end+1} = email_addr;
end

%memory!
command_cell{end+1} = ' --memory ';
command_cell{end+1} = memory;


command = strjoin(command_cell);

%Submit job
[status,result] = system(command);
disp(strtrim(result));


end

function cmd = build_cmd_string(fcn_handle,args)

cmd = [func2str(fcn_handle),'('];
for i=1:numel(args)
    arg = args{i};
    if iscell(arg)
        cmd = [cmd,build_cell_str(arg)];
    elseif ischar(arg)
        cmd = [cmd,sprintf('''%s'',',arg)];
    elseif isscalar(arg)
        cmd = [cmd,sprintf('%g,',arg)];
    elseif isvector(arg)
        cmd = [cmd,build_vec_str(arg)];
    else
        error('Unsupported Argument type, consider saving your arguments to a file and loading them in the function');
    end
end
%Now we should have an extra comma and we need to close the function call
cmd = [cmd(1:end-1),')'];

end

function str_arg = build_cell_str(arg)
str_arg = '{';
for j=1:numel(arg)
    if ischar(arg{j})
        str_arg = [str_arg,sprintf('''%s'',',arg{j})];
    elseif iscell(arg{j})
        str_arg = [str_arg,build_cell_str(arg{j})];
    elseif isscalar(arg{j})
        str_arg = [str_arg,sprintf('%g,',arg{j})];
    elseif isvector(arg{j})
        str_arg = [str_arg,build_vec_str(arg{j})];
    else
        error('Unsupported Argument type, consider saving your arguments to a file and loading them in the function');
    end
end
%Now we should have an extra comma and we need to close the cell
%array
str_arg = [str_arg(1:end-1),'},'];
end

function str_arg = build_vec_str(arg)
str_arg = '[';
for j=1:numel(arg)
    str_arg = [str_arg,sprintf('%g,',arg(j))];
end
%Now we should have an extra comma and we need to close the vector
%array
str_arg = [str_arg(1:end-1),'],'];
end