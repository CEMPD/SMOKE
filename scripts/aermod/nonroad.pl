#!/usr/bin/perl

use strict;
use warnings;

use Text::CSV ();
use Geo::Coordinates::UTM qw(latlon_to_utm latlon_to_utm_force_zone);

require 'aermod.subs';
require 'aermod_np.subs';

# check environment variables
foreach my $envvar (qw(REPORT_JAN REPORT_FEB REPORT_MAR REPORT_APR REPORT_MAY REPORT_JUN 
                       REPORT_JUL REPORT_AUG REPORT_SEP REPORT_OCT REPORT_NOV REPORT_DEC 
                       SOURCE_GROUPS GROUP_PARAMS ATPRO_MONTHLY ATPRO_WEEKLY ATPRO_HOURLY OUTPUT_DIR)) {
  die "Environment variable '$envvar' must be set" unless $ENV{$envvar};
}

# open report files
my @in_fh = ();
foreach my $envvar (qw(REPORT_JAN REPORT_FEB REPORT_MAR REPORT_APR REPORT_MAY REPORT_JUN 
                       REPORT_JUL REPORT_AUG REPORT_SEP REPORT_OCT REPORT_NOV REPORT_DEC)) {
  push @in_fh, open_input($ENV{$envvar});
}

# load source group parameters
print "Reading source groups data...\n";
my $group_fh = open_input($ENV{'GROUP_PARAMS'});

my $csv_parser = Text::CSV->new();
my $header = $csv_parser->getline($group_fh);
$csv_parser->column_names(@$header);

my %group_params;
while (my $row = $csv_parser->getline_hr($group_fh)) {
  $group_params{$row->{'run_group'}} = {
    'release_height' => $row->{'release_height'},
    'sigma_z' => $row->{'sigma_z'}
  };
}
close $group_fh;

# load source group/SCC mapping
$group_fh = open_input($ENV{'SOURCE_GROUPS'});

$csv_parser = Text::CSV->new();
$header = $csv_parser->getline($group_fh);
$csv_parser->column_names(@$header);

my %scc_groups;
while (my $row = $csv_parser->getline_hr($group_fh)) {
  # pad SCCs to 20 characters to match SMOKE 4.0 output
  my $scc = sprintf "%020s", $row->{'SCC'};

  if (exists $scc_groups{$scc}) {
    die "Duplicate SCC $row->{'SCC'} in source groups file";
  }
  
  my $run_group = $row->{'Run Group'};
  unless (exists $group_params{$run_group}) {
    die "Unknown run group name $run_group in source group/SCC mapping file";
  }

  $scc_groups{$scc} = $run_group;
  $group_params{$run_group}{'source_group'} = $row->{'source group'};
}
close $group_fh;

# load temporal profiles
print "Reading temporal profiles...\n";
my $prof_file = $ENV{'ATPRO_MONTHLY'};
my %monthly = read_profiles($prof_file, 12);

$prof_file = $ENV{'ATPRO_WEEKLY'};
my %weekly = read_profiles($prof_file, 7);

$prof_file = $ENV{'ATPRO_HOURLY'};
my %daily = read_profiles($prof_file, 24);

# open output files
print "Creating output files...\n";
my $output_dir = $ENV{'OUTPUT_DIR'};

my $loc_fh = open_output("$output_dir/nonroad_locations.csv");
write_location_header($loc_fh);

my $param_fh = open_output("$output_dir/nonroad_area_params.csv");
write_parameter_header($param_fh);

my $tmp_fh = open_output("$output_dir/nonroad_temporal.csv");
write_temporal_header($tmp_fh);

my $x_fh = open_output("$output_dir/nonroad_emis.csv");
print $x_fh "run_group,region_cd,met_cell,src_id,source_group,smoke_name,ann_value,winter,spring,summer,fall\n";

my %headers;
my @pollutants;
my %sources;
my %emissions;

print "Processing AERMOD sources...\n";

my $month = 0;
for my $fh (@in_fh) {
  $month++;
  print "Reading report for month $month...\n";

  while (my $line = <$fh>) {
    chomp $line;
    next if skip_line($line);
  
    my ($is_header, @data) = parse_report_line($line);
  
    if ($is_header) {
      parse_header(\@data, \%headers, \@pollutants, 'SE Longitude');
      next;
    }
    
    # override assigned temporal profiles
    $data[$headers{'Monthly Prf'}] = '262';
    $data[$headers{'Weekly  Prf'}] = '7';
    $data[$headers{'Mon Diu Prf'}] = '26';
    $data[$headers{'Tue Diu Prf'}] = '26';
    $data[$headers{'Wed Diu Prf'}] = '26';
    $data[$headers{'Thu Diu Prf'}] = '26';
    $data[$headers{'Fri Diu Prf'}] = '26';
    $data[$headers{'Sat Diu Prf'}] = '26';
    $data[$headers{'Sun Diu Prf'}] = '26';
    
    # look up run group based on SCC
    my $scc = $data[$headers{'SCC'}];
    unless (exists $scc_groups{$scc}) {
      die "No run group defined for SCC $scc";
    }
    my $run_group = $scc_groups{$scc};
  
    # build cell identifier
    my $cell = "G" . sprintf("%03d", $data[$headers{'X cell'}]) .
               "R" . sprintf("%03d", $data[$headers{'Y cell'}]);

    my $source_id = join(":::", $run_group, $cell);
    unless (exists $sources{$source_id}) {
      $sources{$source_id} = 1;

      my @common;
      push @common, $run_group;
      push @common, $cell;
      push @common, "12_";

      # prepare location output
      my @output = @common;
      my $sw_lat = $data[$headers{'SW Latitude'}];
      my $sw_lon = $data[$headers{'SW Longitude'}];
      my ($zone, $utm_x, $utm_y) = latlon_to_utm(23, $sw_lat, $sw_lon);
      my $outzone = $zone;
      $outzone =~ s/\D//g; # strip latitude band designation from UTM zone
      push @output, "AREAPOLY";
      push @output, $utm_x;
      push @output, $utm_y;
      push @output, $outzone;
      push @output, $sw_lon;
      push @output, $sw_lat;
      print $loc_fh join(',', @output) . "\n";
  
      # prepare parameters output
      @output = @common;
      push @output, $group_params{$run_group}{'release_height'};
      push @output, "4"; # number of vertices
      push @output, $group_params{$run_group}{'sigma_z'};
      push @output, $utm_x;
      push @output, $utm_y;

      my $nw_lat = $data[$headers{'NW Latitude'}];
      my $nw_lon = $data[$headers{'NW Longitude'}];
      ($zone, $utm_x, $utm_y) = latlon_to_utm_force_zone(23, $zone, $nw_lat, $nw_lon);
      push @output, $utm_x;
      push @output, $utm_y;

      my $ne_lat = $data[$headers{'NE Latitude'}];
      my $ne_lon = $data[$headers{'NE Longitude'}];
      ($zone, $utm_x, $utm_y) = latlon_to_utm_force_zone(23, $zone, $ne_lat, $ne_lon);
      push @output, $utm_x;
      push @output, $utm_y;

      my $se_lat = $data[$headers{'SE Latitude'}];
      my $se_lon = $data[$headers{'SE Longitude'}];
      ($zone, $utm_x, $utm_y) = latlon_to_utm_force_zone(23, $zone, $se_lat, $se_lon);
      push @output, $utm_x;
      push @output, $utm_y;

      push @output, $sw_lon;
      push @output, $sw_lat;
      push @output, $nw_lon;
      push @output, $nw_lat;
      push @output, $ne_lon;
      push @output, $ne_lat;
      push @output, $se_lon;
      push @output, $se_lat;
      print $param_fh join(',', @output) . "\n";
    
      # prepare temporal profile output
      my ($qflag, @factors) = get_factors(\%headers, \@data, \%monthly, \%weekly, \%daily);
    
      @output = @common;
      push @output, $qflag;
      push @output, map { sprintf('%.4f', $_) } @factors;
      print $tmp_fh join(',', @output) . "\n";
    }
  
    # store emissions by season
    my $season = 'winter';  # dec, jan, feb
    if ($month >= 3 && $month <= 5) {  # mar, apr, may
      $season = 'spring';
    } elsif ($month >= 6 && $month <= 8) {  # jun, jul, aug
      $season = 'summer';
    } elsif ($month >= 9 && $month <= 11) {  # sep, oct, nov
      $season = 'fall';
    }
    
    my $region = $data[$headers{'Region'}];
    $region = substr($region, -5);
    foreach my $poll (@pollutants) {
      my $emis_id = join(":::", $source_id, $region, $poll);
      unless (exists $emissions{$emis_id}) {
        $emissions{$emis_id} = {
          'winter' => 0,
          'spring' => 0,
          'summer' => 0,
          'fall' => 0
        };
      }
      $emissions{$emis_id}{$season} = 
        $emissions{$emis_id}{$season} + $data[$headers{$poll}];
    }
  }
}

# prepare crosswalk output
for my $emis_id (keys %emissions) {
  my ($run_group, $cell, $region, $poll) = split(/:::/, $emis_id);

  my @output;
  push @output, $run_group;
  push @output, $region;
  push @output, $cell;
  push @output, "12_1";
  push @output, '"'.$group_params{$run_group}{'source_group'}.'"';
  push @output, $poll;
  
  my $total = $emissions{$emis_id}{'winter'} +
              $emissions{$emis_id}{'spring'} +
              $emissions{$emis_id}{'summer'} +
              $emissions{$emis_id}{'fall'};
  push @output, $total;
  push @output, $emissions{$emis_id}{'winter'};
  push @output, $emissions{$emis_id}{'spring'};
  push @output, $emissions{$emis_id}{'summer'};
  push @output, $emissions{$emis_id}{'fall'};
  print $x_fh join(',', @output) . "\n";
}

for my $fh (@in_fh) {
  close $fh;
}
close $loc_fh;
close $param_fh;
close $tmp_fh;
close $x_fh;

print "Done.\n";
