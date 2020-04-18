#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use Cwd qw[abs_path];


my $scriptpath = abs_path($0);

sub usage_exit
{
	print("Usage: $0 [options]\n");
	print("Run the a set of jobs on the cluster.\n");
	print("  -h, --help         Displays this information\n");
	print("  -b, --batch        Batch Filename\n");
	print("  -a, --array        Array range\n");
	print("  -d, --modules       modules tool\n");
	print("  -p, --partition        Parition name: defaults to 'campus' always.\n");
	print("  -m, --memory				Memory in GB (e.g. 4 for 4000MB or 4GB)\n");
	print("  -r, --runtime			Runtime in hours; default 2\n");
	print("  -c, --cores				Cores (e.g. 2)\n");
	print("  -l, --log          Log files\n");
	print("  -o, --slurm_options      Additional slurm options. Optional\n");
	print("  -s, --script           Bash script with slurm options\n");
	exit;
}

my $help;
my $submit_host = `hostname`;
chomp($submit_host);
my $batch_filename;
my $parition_name = "campus";
my $modules = "R/3.6.0-foss-2016b-fh1";
my $array = "";
my $launch = 0;
my $memory = 4;
my $runtime = "1";
my $cores = 1;
my $nodes = 1;
my $slurm_options;
my $script = "";
my $logDir = "slurm_logs/";

GetOptions
(
    'help'         				=> \$help,
    'batch|b=s'      				=> \$batch_filename,
    'parition|p=s'    			=> \$parition_name,
    'slurm_options|o'			=> \$slurm_options,
    'log|l=s'					=> \$logDir,
    'memory|m=s' 						=> \$memory,
    'cores|c=i'							=> \$cores,
    'nodes=i'                 => \$nodes,
    'runtime|r=s'						=> \$runtime,
    'loadModules|d=s'						=> \$modules,
    'array|a=s'							=> \$array,
    'script=s'                  => \$script
);

not defined $help or usage_exit();
defined $batch_filename or usage_exit();

$runtime = "${runtime}:0:0";
$memory = "${memory}000";
mkdir $logDir;
if ($modules ne ""){
    $modules = "ml " . $modules;
}

if ($script eq ""){
	($script = $batch_filename) =~ s/.sh|.txt/.slurm/g;
}

print "PARAMETERS:\nbatch=$batch_filename\nsubmission_host=$submit_host\nmemory=$memory\ncores=$cores\nnodes=$nodes\nruntime=$runtime\nlogDir=$logDir\nslurm_options=$slurm_options\nparition_name=$parition_name\nmodules=$modules\narray=$array\nlaunch=$launch\n";

$batch_filename = abs_path($batch_filename);
my $batch_jobname_prefix = basename($batch_filename);
#(my $batch_jobname_prefix = basename($batch_filename)) =~ s/.sh//g;

## output bash script with slurm parameters/options
open SCRIPT, ">$script" || die "Can't write to $script.\n";

print SCRIPT "#!/bin/bash\n\n";
print SCRIPT "#SBATCH --job-name=${batch_jobname_prefix}\n";
print SCRIPT "#SBATCH --partition=${parition_name}\n";
print SCRIPT "#SBATCH --nodes=${nodes}\n";
print SCRIPT "#SBATCH --cpus-per-task=${cores}\n";
print SCRIPT "#SBATCH --mem=${memory}\n";
print SCRIPT "#SBATCH --time=${runtime}\n";
print SCRIPT "#SBATCH -o ${logDir}/${batch_jobname_prefix}.o.%A.%a\n";
print SCRIPT "#SBATCH $slurm_options\n" if ($slurm_options);
print SCRIPT "\n";

my $module_env = "export PATH=/app/bin:\$PATH
source /app/Lmod/lmod/lmod/init/bash
export LMOD_PACKAGE_PATH=/app/Lmod
module use /app/easybuild/modules/all
module use /app/Modules/modulefiles";
print SCRIPT "$module_env\n$modules\n\n";

#print SCRIPT "cmd=\`sed -n \"\${SLURM_ARRAY_TASK_ID}{p;q;}\" $batch_filename\`\n";
print SCRIPT "eval \`sed -n \"\${SLURM_ARRAY_TASK_ID}{p;q;}\" $batch_filename\`\n";

close SCRIPT;

#my $jobs_to_start = scalar @jobs;
#print "Starting $jobs_to_start jobs\n";
if ($array eq ""){
	$array = `cat $batch_filename | wc -l`; chomp($array);
	$array = "1-$array";
}

my $sbatch_job = "sbatch --array=$array $script\n";
print "\n" . $sbatch_job . "\n";


