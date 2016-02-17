#! /usr/bin/perl -w
#
# Fortran 90/95 dependency checker.

use 5.006_000;
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Getopt::Long;
use Cwd;

sub modules_used_in {
  my $file = shift;
  my @excluded = @_;
  local(*FILE);
  local($_);

  if (-f $file) {
    open(FILE, $file) || warn "Can't open $file: $!\n";
  } else {
    print STDERR "Warning: file \"${file}\" does not exist\".\n";
    return;
  }

  my @modules = ();
  my $i = 0;

  while (<FILE>) {
    # Fortran 90 "use"
    if (/^\s*use\s+(\w+)/i) {   
      my $modname = lc($1);
      next if file_matches_module_name($file,$modname);
      next if mod_not_included($modname,@excluded);
      $modules[$i++] = $modname;
    }
  }
  uniq(@modules);
}

sub file_matches_module_name {
    my $file = shift;
    my $modname = shift;

    ($file =~ m/^(.*\/|)$modname\.[FfCc][9Uu][0Ff]$/i?1:0);
}

sub uniq {
  my @in = @_;
  undef my %saw;
  return(grep(!$saw{$_}++, @in));
}

sub mod_not_included {
  my $mod = shift;
  my @excluded_mods = @_;
  my $nc;
  my @not_include = qw(f90_unix nas_system iso_c_binding
                       synchronize_variables 
                       f90_unix_io
                       f90_unix_dir 
                       f90_unix_proc
                       punditinterface
                       smemrd_common
                       smemrd_interface
                       zoltan);
  push @not_include, @excluded_mods;
  foreach $nc (@not_include) {
    if ( $mod eq $nc ) {
      return 1
    }
  }
}

sub files_included_in {
  my $file = shift;
  local(*FILE);
  local($_);

  if (-f $file) {
    open(FILE, $file) || warn "Can't open $file: $!\n";
  } else {
    return;
  }

  my @files = ();
  my $i = 0;

#print STDERR "Searching for includes in \"$file\".\n";
  while (<FILE>) {
    # Fortran 90 "include"
#print STDERR "\"$file\" line is \"$_\".\n";
    if (/^\s*include\s+\'(.+)\'/i) {   
#print STDERR "include in \"$file\" as \"$_\".\n";
      $files[$i++] = main::basename($1);
    } else {
#print STDOUT "$_";
    }
  }
  uniq(@files);
}

sub modules_defined_in {
  my $file = shift;
  local(*FILE);
  local($_);

  if (-f $file) {
    open(FILE, $file) || warn "Can't open $file: $!\n";
  } else {
    return;
  }

  my @modules = ();
  my $i = 0;

  while (<FILE>) {
    # Fortran 90 module
    if (/^\s*module\s+(\w+)/i) {
      my $modname = lc($1);
      next if (lc($modname) eq "procedure");
      $modules[$i++] = lc($1);
    }
  }
  uniq(@modules);
}

sub fortran_files_within {
  my @searchPath = @_;
  my @source = map(glob("$_/*.[Ff][9Uu][0Ff]"),@searchPath);
}

sub file_search {
  my $file = $_[0];
  my $dir = $_[1];
  my $file_found;

  opendir (DIR, "$dir") or die "Failed to open directory: $dir\n";
  $file_found = grep { /^$file$/ } readdir DIR;
  closedir (DIR);
  #print "dir ::@dirs_found\n";
  return $file_found;
}

sub strip_extensions {
  my @files = @_;
  my @modules = map { (my $file = $_) =~ s/\..*$//; $file } @files;
  return @modules;
}

my ($mod_ext, @inlined_dirs, @inc_dirs, @excluded_files);
$mod_ext = "mod";

Getopt::Long::Configure( "bundling" );
GetOptions("I=s@"      => \@inc_dirs, 
           "mod_ext=s" => \$mod_ext,   "L=s@"  => \@inlined_dirs,
           "X=s@" => \@excluded_files);

my $source = $ARGV[0];
my $srcdir = $ARGV[1];

my @excluded_mods = strip_extensions(@excluded_files);
my @searchPath = @inc_dirs;

my @modules = modules_used_in($source,@excluded_mods);
my @included_files = files_included_in($source);
my @mod_name = modules_defined_in($source);

#@modules = remove_mods_defined_in_source(\@mod_name,\@modules);

my ($file, @included_files_path);

#search through inlined files for modules
foreach $file (@included_files) {
  foreach my $path (@inlined_dirs) {
    if (file_search($file,$path)) {
      push(@included_files_path,"$path/$file");
      push(@modules,modules_used_in("$path/$file",@excluded_mods));
      last;
    };
  };
};

@modules = uniq(@modules);

# Remove overlap if modules are defined in the same file
# Fixed issues with some files in Utils subdir

my @extra_mods;

foreach my $m (@mod_name) {
  push (@extra_mods, $m) unless file_matches_module_name($source,$m);
}

my @diff;
my %repeats;
  
for (@modules, @extra_mods) { $repeats{$_}++ }
for (keys %repeats) {
  push @diff, $_ unless $repeats{$_} > 1;
}

@modules = @diff;

my %module_dependency = map{$_ => "module not found"} @modules;

#get the currect directory - will check to see if the directory is a complex dir later
my $cwd = basename(getcwd());

foreach my $mod (keys %module_dependency) {
  if(file_search("$mod.f90",$srcdir) || file_search("$mod.F90",$srcdir) ||
     file_search("$mod.cuf",$srcdir) || file_search("$mod.CUF",$srcdir) || "Complex" eq $cwd ) {
    $module_dependency{$mod} = "$mod.$mod_ext"
  };
};

foreach my $path (@searchPath) {
  foreach my $mod (keys %module_dependency) {
    if(file_search("$mod.$mod_ext",$path)) {
      $module_dependency{$mod} = "$path/$mod.$mod_ext"
    };
  };
};

my ($dot_o, @dot_mod);

($dot_o = $source) =~ s/\.[Ff]90/\.o/g;
$dot_o =~ s/\.[Cc][Uu][Ff]/\.o/g;
$dot_o =~ s/^.*\///g;
@dot_mod = map("$_.$mod_ext",@mod_name);

(my $f90_file = $source) =~ s/^.*\///g;

my @output;
push(@output, $dot_o);
push(@output, @dot_mod);
push(@output, ":");
push(@output, $f90_file);
push(@output, @included_files_path);

foreach my $mod (keys %module_dependency) {
  if ($module_dependency{$mod} eq "module not found") {
    print STDERR "Warning: unable to find module $mod for file $source\n";
  } else {
    push(@output,$module_dependency{$mod});
  }
}; 

print STDOUT "@output\n";

foreach my $mod (keys %module_dependency) {
  print STDOUT "$module_dependency{$mod} :\n";
}; 
foreach my $include (@included_files_path) {
  print STDOUT "$include :\n"
};
