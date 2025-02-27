#!/usr/bin/perl

use strict;
use warnings;

use Text::CSV ();
use Geo::Coordinates::UTM qw(latlon_to_utm latlon_to_utm_force_zone);
use Date::Simple qw(ymd);

require 'aermod.subs';
require 'aermod_np.subs';

my $grid_prefix = $ENV{'GRID_PREFIX'} || '12_';
my $run_group_suffix = $ENV{'RUN_GROUP_SUFFIX'} || '12';

my @days_in_month = (0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31);

my $rwc_run_group;

# check environment variables
foreach my $envvar (qw(REPORT SOURCE_GROUPS GROUP_PARAMS YEAR
                       ATPRO_MONTHLY ATPRO_DAILY ATPRO_WEEKLY ATPRO_HOURLY OUTPUT_DIR)) {
  die "Environment variable '$envvar' must be set" unless $ENV{$envvar};
}

my $year = $ENV{'YEAR'};

# open report file
my $in_fh = open_input($ENV{'REPORT'});

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
  my $scc = sprintf "%020s", $row->{'scc'};

  if (exists $scc_groups{$scc}) {
    die "Duplicate SCC $row->{'scc'} in source groups file";
  }
  
  my $run_group = $row->{'run_group'};
  unless (exists $group_params{$run_group}) {
    die "Unknown run group name $run_group in source group/SCC mapping file";
  }

  $scc_groups{$scc}{'run_group'} = $run_group;
  $scc_groups{$scc}{'source_group'} = $row->{'source_group'};
}
close $group_fh;

# load temporal profiles
print "Reading temporal profiles...\n";
my $prof_file = $ENV{'ATPRO_MONTHLY'};
my %monthly = read_profiles($prof_file, 12);

$prof_file = $ENV{'ATPRO_DAILY'};
my %dayofmonth = read_dom_profiles($prof_file, 31, \@days_in_month);

$prof_file = $ENV{'ATPRO_WEEKLY'};
my %weekly = read_profiles($prof_file, 7);

$prof_file = $ENV{'ATPRO_HOURLY'};
my %daily = read_profiles($prof_file, 24);

# check output directories
print "Checking output directories...\n";
my $output_dir = $ENV{'OUTPUT_DIR'};
die "Missing output directory $output_dir" unless -d $output_dir;

foreach my $dir (qw(locations parameters temporal emis xwalk)) {
  die "Missing output directory $output_dir/$dir" unless -d $output_dir . '/' . $dir;
}

my %handles;

my %headers;
my @pollutants;
my %sources;
my %gridded_emissions;  # emissions by grid cell, source group, county, and pollutant

# RWC-specific collections
my %county_emissions;  # sum of all emissions by county and SCC
my %hour_of_year_factors;  # hour of year factors by county and SCC

print "Processing AERMOD sources...\n";
while (my $line = <$in_fh>) {
  chomp $line;
  next if skip_line($line);

  my ($is_header, @data) = parse_report_line($line);

  if ($is_header) {
    parse_header(\@data, \%headers, \@pollutants, 'SE Longitude');
    next;
  }

  # look up run group based on SCC
  my $scc = $data[$headers{'SCC'}];
  unless (exists $scc_groups{$scc}) {
    die "No run group defined for SCC $scc";
  }
  
  # all SCCs must use the same run group, so save the first one seen and check subsequent SCCs
  $rwc_run_group = $scc_groups{$scc}{'run_group'} unless $rwc_run_group;
  unless ($scc_groups{$scc}{'run_group'} eq $rwc_run_group) {
    die "SCC $scc uses different run group '$scc_groups{$scc}{'run_group'}' than default group '$rwc_run_group'";
  }
  
  my $source_group = $scc_groups{$scc}{'source_group'};

  # build cell identifier
  my $cell = "G" . sprintf("%03d", $data[$headers{'X cell'}]) .
             "R" . sprintf("%03d", $data[$headers{'Y cell'}]);

  unless (exists $sources{$cell}) {
    $sources{$cell} = 1;

    my @common;
    push @common, $rwc_run_group . $run_group_suffix;
    push @common, $cell;
    push @common, "${grid_prefix}1";

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
    my $file = "$output_dir/locations/${rwc_run_group}${run_group_suffix}_locations.csv";
    unless (exists $handles{$file}) {
      my $fh = open_output($file);
      write_location_header($fh);
      $handles{$file} = $fh;
    }
    my $loc_fh = $handles{$file};
    print $loc_fh join(',', @output) . "\n";

    # prepare parameters output
    @output = @common;
    push @output, $group_params{$rwc_run_group}{'release_height'};
    push @output, "4"; # number of vertices
    push @output, $group_params{$rwc_run_group}{'sigma_z'};
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
    $file = "$output_dir/parameters/${rwc_run_group}${run_group_suffix}_area_params.csv";
    unless (exists $handles{$file}) {
      my $fh = open_output($file);
      write_parameter_header($fh);
      $handles{$file} = $fh;
    }
    my $param_fh = $handles{$file};
    print $param_fh join(',', @output) . "\n";
  }

  # store emissions
  my $region = $data[$headers{'Region'}];
  $region = substr($region, -5);
  
  unless (exists $county_emissions{$region}{'all'}) {
    $county_emissions{$region}{'all'} = 0;
  }
  unless (exists $county_emissions{$region}{$scc}) {
    $county_emissions{$region}{$scc} = 0;
  }
  unless (exists $gridded_emissions{$cell}{$source_group}{$region}{'all'}) {
    $gridded_emissions{$cell}{$source_group}{$region}{'all'} = 0;
  }
  
  foreach my $poll (@pollutants) {
    unless (exists $gridded_emissions{$cell}{$source_group}{$region}{$poll}) {
      $gridded_emissions{$cell}{$source_group}{$region}{$poll} = 0;
    }
    my $emis_value = $data[$headers{$poll}];
    
    $gridded_emissions{$cell}{$source_group}{$region}{$poll} += $emis_value;
    $gridded_emissions{$cell}{$source_group}{$region}{'all'} += $emis_value;
    
    $county_emissions{$region}{$scc} += $emis_value;
    $county_emissions{$region}{'all'} += $emis_value;
  }
  
  # calculate hour of year factors for county and SCC
  unless (exists $hour_of_year_factors{$region}{$scc}) {
    my $monthly_prof = $data[$headers{'Monthly Prf'}];
    die "Unknown monthly profile code: $monthly_prof" unless exists $monthly{$monthly_prof};
    my @monthly_factors = @{$monthly{$monthly_prof}};
    
    my $dom_prof = $data[$headers{'Day-Month Prf'}];
    my %dom_factors;
    if ($dom_prof) {
      die "Unknown day-of-month profile code: $dom_prof" unless exists $dayofmonth{$dom_prof};
      %dom_factors = %{$dayofmonth{$dom_prof}};
    }
    
    my $weekly_prof = $data[$headers{'Weekly  Prf'}];
    my @weekly_factors;
    if ($weekly_prof) {
      die "Unknown weekly profile code: $weekly_prof" unless exists $weekly{$weekly_prof};
      @weekly_factors = @{$weekly{$weekly_prof}};
    }
    
    my $monday_prof = $data[$headers{'Mon Diu Prf'}];
    # check that all days of the week use the same profile
    unless ($monday_prof eq $data[$headers{'Tue Diu Prf'}] &&
            $monday_prof eq $data[$headers{'Wed Diu Prf'}] &&
            $monday_prof eq $data[$headers{'Thu Diu Prf'}] &&
            $monday_prof eq $data[$headers{'Fri Diu Prf'}] &&
            $monday_prof eq $data[$headers{'Sat Diu Prf'}] &&
            $monday_prof eq $data[$headers{'Sun Diu Prf'}]) {
      die "All days of the week must use the same hourly temporal profile";
    }
    die "Unknown hourly profile code: $monday_prof" unless exists $daily{$monday_prof};
    my @hourly_factors = @{$daily{$monday_prof}};
    
    my @factors;
    my $month = 1;
    # determine starting day of the week (convert from 0 = Sunday to 0 = Monday for temporal profiles)
    my $day_of_week = (ymd($year, $month, 1)->day_of_week - 1) % 7;
    foreach my $month_factor (@monthly_factors) {
      # note: factors average to 1 rather than summing to 1, so adjust when fractions are needed
      if ($dom_prof) {
        my $day_of_month = 1;
        foreach my $dom_factor (@{$dom_factors{$month}}) {
          last if $day_of_month > $days_in_month[$month];
          push @factors, map { $_/24 * $dom_factor/$days_in_month[$month] * $month_factor/12 } @hourly_factors;
          $day_of_month++;
        }
      } else {
        my $month_to_week = 7 / $days_in_month[$month];
        for (my $day_of_month = 1; $day_of_month <= $days_in_month[$month]; $day_of_month++) {
          push @factors, map { $_/24 * $weekly_factors[$day_of_week]/7 * $month_to_week * $month_factor/12 } @hourly_factors;
          $day_of_week = ($day_of_week + 1) % 7;
        }
      }
      $month++;
    }
    
    $hour_of_year_factors{$region}{$scc} = \@factors;
  }
}

# prepare temporal profile output
for my $region (sort keys %county_emissions) {
  my $state = substr($region, 0, 2);
  my $file = "$output_dir/temporal/${rwc_run_group}${run_group_suffix}_${state}_hourly.csv";
  unless (exists $handles{$file}) {
    my $fh = open_output($file);
    print $fh "run_group,region_cd,year,month,day,hour,factor\n";
    $handles{$file} = $fh;
  }
  my $tmp_fh = $handles{$file};

  my $total_emissions = $county_emissions{$region}{'all'};
  
  my @common;
  push @common, $rwc_run_group . $run_group_suffix;
  push @common, $region;
  push @common, substr($year, -2);
  
  my $index = 0;
  for (my $month = 1; $month <= 12; $month++) {
    for (my $day_of_month = 1; $day_of_month <= $days_in_month[$month]; $day_of_month++) {
      for (my $hour = 1; $hour <= 24; $hour++) {
        my $hourly_emissions = 0;
        for my $scc (keys %{$county_emissions{$region}}) {
          next if $scc eq 'all';
          
          my $factors_ref = $hour_of_year_factors{$region}{$scc};
          my $size = scalar(@$factors_ref);
          die "$region, $scc, $size" unless $size == 8760;
          $hourly_emissions += $county_emissions{$region}{$scc} * $factors_ref->[$index];
        }
      
        my @output = @common;
        push @output, $month;
        push @output, $day_of_month;
        push @output, $hour;
        push @output, sprintf('%.8f', 8760 * $hourly_emissions / $total_emissions);
        print $tmp_fh join(',', @output) . "\n";
        
        $index++;
      }
    }
  }
}
  
# prepare crosswalk output and grid cell / county assignment
my $cnty_fh = open_output("$output_dir/xwalk/${rwc_run_group}${run_group_suffix}_county-to-gridcell.csv");
print $cnty_fh "run_group,region_cd,met_cell,src_id\n";

my $x_fh = open_output("$output_dir/emis/${grid_prefix}${rwc_run_group}${run_group_suffix}_emis.csv");
print $x_fh "run_group,region_cd,met_cell,src_id,source_group,smoke_name,ann_value\n";

for my $cell (sort keys %gridded_emissions) {
  my $cell_assignment;
  my $prev_emis = -999;
  for my $source_group (sort keys %{$gridded_emissions{$cell}}) {
    for my $region (sort keys %{$gridded_emissions{$cell}{$source_group}}) {
      for my $poll (sort keys %{$gridded_emissions{$cell}{$source_group}{$region}}) {
        next if $poll eq 'all';
        next if $gridded_emissions{$cell}{$source_group}{$region}{$poll} == 0.0;

        my @output;
        push @output, $rwc_run_group . $run_group_suffix;
        push @output, $region;
        push @output, $cell;
        push @output, "${grid_prefix}1";
        push @output, $source_group;
        push @output, $poll;
        push @output, $gridded_emissions{$cell}{$source_group}{$region}{$poll};
        print $x_fh join(',', @output) . "\n";
      }
    
      # check if current county has greater emissions than previous for this grid cell
      if ($gridded_emissions{$cell}{$source_group}{$region}{'all'} > $prev_emis) {
        $cell_assignment = $region;
        $prev_emis = $gridded_emissions{$cell}{$source_group}{$region}{'all'};
      }
    }
  }
  
  my @output;
  push @output, $rwc_run_group . $run_group_suffix;
  push @output, $cell_assignment;
  push @output, $cell;
  push @output, "${grid_prefix}1";
  print $cnty_fh join(',', @output) . "\n";
}

close $in_fh;
foreach my $fh (values %handles) {
  close $fh;
}
close $cnty_fh;
close $x_fh;

print "Done.\n";
