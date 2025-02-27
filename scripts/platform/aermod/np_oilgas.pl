#!/usr/bin/perl

use strict;
use warnings;

use Text::CSV ();
use Geo::Coordinates::UTM qw(latlon_to_utm latlon_to_utm_force_zone);

require 'aermod.subs';
require 'aermod_np.subs';

my $grid_prefix = $ENV{'GRID_PREFIX'} || '4_';
my $run_group_suffix = $ENV{'RUN_GROUP_SUFFIX'} || '4';

# set up grid cell subsetting
my $num_cells = 1;
unless ($ENV{'USE_GRID_CELL_SUBSETTING_YN'} &&
        $ENV{'USE_GRID_CELL_SUBSETTING_YN'} eq 'N') {
  my $outer_size = $ENV{'OUTER_GRID_CELL_SIZE'} || 12;
  my $inner_size = $ENV{'INNER_GRID_CELL_SIZE'} || 4;
  
  # check that outer size is a multiple of inner size
  if ($outer_size % $inner_size != 0) {
    die "OUTER_GRID_CELL_SIZE must be a multiple of INNER_GRID_CELL_SIZE";
  }
  
  $num_cells = $outer_size / $inner_size;
};

my $oilgas_run_group;

# check environment variables
foreach my $envvar (qw(REPORT SOURCE_GROUPS GROUP_PARAMS 
                       ATPRO_MONTHLY ATPRO_WEEKLY ATPRO_HOURLY OUTPUT_DIR)) {
  die "Environment variable '$envvar' must be set" unless $ENV{$envvar};
}

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

$prof_file = $ENV{'ATPRO_WEEKLY'};
my %weekly = read_profiles($prof_file, 7);

$prof_file = $ENV{'ATPRO_HOURLY'};
my %daily = read_profiles($prof_file, 24);

# check output directories
print "Checking output directories...\n";
my $output_dir = $ENV{'OUTPUT_DIR'};
die "Missing output directory $output_dir" unless -d $output_dir;

foreach my $dir (qw(locations parameters temporal emis)) {
  die "Missing output directory $output_dir/$dir" unless -d $output_dir . '/' . $dir;
}

my %handles;

my %headers;
my @pollutants;
my %sources;
my %gridded_emissions;  # emissions by grid cell, source group, county, and pollutant
my %county_emissions;  # sum of all emissions by county and SCC
my %month_of_year_factors;  # month of year factors by county and SCC

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
  $oilgas_run_group = $scc_groups{$scc}{'run_group'} unless $oilgas_run_group;
  unless ($scc_groups{$scc}{'run_group'} eq $oilgas_run_group) {
    die "SCC $scc uses different run group '$scc_groups{$scc}{'run_group'}' than default group '$oilgas_run_group'";
  }
  
  my $source_group = $scc_groups{$scc}{'source_group'};

  # build cell identifier
  my $resolution_id = '1';
  my $cell;
  if ($num_cells > 1) {
    my $x_inner = $data[$headers{'X cell'}];
    my $y_inner = $data[$headers{'Y cell'}];
    my $x_outer = int(($x_inner + $num_cells - 1) / $num_cells);
    my $y_outer = int(($y_inner + $num_cells - 1) / $num_cells);
    $cell = "G" . sprintf("%03d", $x_outer) .
            "R" . sprintf("%03d", $y_outer);
  
    $resolution_id = ((($x_inner - 1) % $num_cells) + 1) +
                     ((($y_inner - 1) % $num_cells) * $num_cells);
  } else {
    $cell = "G" . sprintf("%03d", $data[$headers{'X cell'}]) .
            "R" . sprintf("%03d", $data[$headers{'Y cell'}]);
  }
  
  my $source_id = join(":::", $cell, $resolution_id);
  unless (exists $sources{$source_id}) {
    $sources{$source_id} = 1;

    my @common;
    push @common, $oilgas_run_group . $run_group_suffix;
    push @common, $cell;
    push @common, "${grid_prefix}${resolution_id}";

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
    my $file = "$output_dir/locations/${oilgas_run_group}${run_group_suffix}_locations.csv";
    unless (exists $handles{$file}) {
      my $fh = open_output($file);
      write_location_header($fh);
      $handles{$file} = $fh;
    }
    my $loc_fh = $handles{$file};
    print $loc_fh join(',', @output) . "\n";

    # prepare parameters output
    @output = @common;
    push @output, $group_params{$oilgas_run_group}{'release_height'};
    push @output, "4"; # number of vertices
    push @output, $group_params{$oilgas_run_group}{'sigma_z'};
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
    $file = "$output_dir/parameters/${oilgas_run_group}${run_group_suffix}_area_params.csv";
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
  unless (exists $gridded_emissions{$source_id}{$source_group}{$region}{'all'}) {
    $gridded_emissions{$source_id}{$source_group}{$region}{'all'} = 0;
  }
  
  foreach my $poll (@pollutants) {
    unless (exists $gridded_emissions{$source_id}{$source_group}{$region}{$poll}) {
      $gridded_emissions{$source_id}{$source_group}{$region}{$poll} = 0;
    }
    my $emis_value = $data[$headers{$poll}];
    
    $gridded_emissions{$source_id}{$source_group}{$region}{$poll} += $emis_value;
    $gridded_emissions{$source_id}{$source_group}{$region}{'all'} += $emis_value;
    
    $county_emissions{$region}{$scc} += $emis_value;
    $county_emissions{$region}{'all'} += $emis_value;
  }
  
  # store month of year factors for county and SCC
  unless (exists $month_of_year_factors{$region}{$scc}) {
    my ($qflag, @factors) = get_factors(\%headers, \@data, \%monthly, \%weekly, \%daily);
    die "SCC $scc in region $region uses non-monthly factors" if ($qflag ne 'MONTH');
    $month_of_year_factors{$region}{$scc} = \@factors;
  }
}

# calculate county-wide temporal profiles
for my $region (sort keys %county_emissions) {
  $month_of_year_factors{$region}{'all'} = [];
  my $total_emissions = $county_emissions{$region}{'all'};
  
  for (my $month = 0; $month < 12; $month++) {
    my $monthly_emissions = 0;
    for my $scc (keys %{$county_emissions{$region}}) {
      next if $scc eq 'all';
    
      my $factors_ref = $month_of_year_factors{$region}{$scc};
      $monthly_emissions += $county_emissions{$region}{$scc} * $factors_ref->[$month] / 12;
    }
    
    $month_of_year_factors{$region}{'all'}[$month] = 12 * $monthly_emissions / $total_emissions;
  }
}

# prepare temporal profile and crosswalk output
my $tmp_fh = open_output("$output_dir/temporal/${oilgas_run_group}${run_group_suffix}_temporal.csv");
write_temporal_header($tmp_fh);

my $x_fh = open_output("$output_dir/emis/${grid_prefix}${oilgas_run_group}${run_group_suffix}_emis.csv");
print $x_fh "run_group,region_cd,met_cell,src_id,source_group,smoke_name,ann_value\n";

for my $source_id (sort keys %gridded_emissions) {
  my ($cell, $resolution_id) = split(/:::/, $source_id);
  my $cell_assignment;
  my $prev_emis = -999;
  for my $source_group (sort keys %{$gridded_emissions{$source_id}}) {
    for my $region (sort keys %{$gridded_emissions{$source_id}{$source_group}}) {
      for my $poll (sort keys %{$gridded_emissions{$source_id}{$source_group}{$region}}) {
        next if $poll eq 'all';
        next if $gridded_emissions{$source_id}{$source_group}{$region}{$poll} == 0.0;

        my @output;
        push @output, $oilgas_run_group . $run_group_suffix;
        push @output, $region;
        push @output, $cell;
        push @output, "${grid_prefix}${resolution_id}";
        push @output, $source_group;
        push @output, $poll;
        push @output, $gridded_emissions{$source_id}{$source_group}{$region}{$poll};
        print $x_fh join(',', @output) . "\n";
      }
      
      # check if current county has greater emissions than previous for this grid cell
      if ($gridded_emissions{$source_id}{$source_group}{$region}{'all'} > $prev_emis) {
        $cell_assignment = $region;
        $prev_emis = $gridded_emissions{$source_id}{$source_group}{$region}{'all'};
      }
    }
  }
  
  my @output;
  push @output, $oilgas_run_group . $run_group_suffix;
  push @output, $cell;
  push @output, "${grid_prefix}${resolution_id}";
  push @output, 'MONTH';
  push @output, map { sprintf('%.8f', $_) } @{ $month_of_year_factors{$cell_assignment}{'all'} };
  print $tmp_fh join(',', @output) . "\n";
}

close $in_fh;
foreach my $fh (values %handles) {
  close $fh;
}
close $tmp_fh;
close $x_fh;

print "Done.\n";
