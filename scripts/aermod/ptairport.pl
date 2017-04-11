#!/usr/bin/perl

use strict;
use warnings;

use Text::CSV ();
use Geo::Coordinates::UTM qw(latlon_to_utm);

require 'aermod.subs';
require 'aermod_pt.subs';

# check environment variables
foreach my $envvar (qw(REPORT RUNWAYS PTPRO_MONTHLY PTPRO_WEEKLY PTPRO_HOURLY OUTPUT_DIR)) {
  die "Environment variable '$envvar' must be set" unless $ENV{$envvar};
}

# open report file
my $in_fh = open_input($ENV{'REPORT'});

# load runway data
print "Reading runway data...\n";
my $runway_fh = open_input($ENV{'RUNWAYS'});

my $csv_parser = Text::CSV->new();
my $header = $csv_parser->getline($runway_fh);
$csv_parser->column_names(@$header);

my %runways;
while (my $row = $csv_parser->getline_hr($runway_fh)) {
  unless (exists $runways{$row->{'facility_id'}}) {
    $runways{$row->{'facility_id'}} = [];
  }
  push @{$runways{$row->{'facility_id'}}}, $row;
}
close $runway_fh;

# load temporal profiles
print "Reading temporal profiles...\n";
my $prof_file = $ENV{'PTPRO_MONTHLY'};
my %monthly = read_profiles($prof_file, 12);

$prof_file = $ENV{'PTPRO_WEEKLY'};
my %weekly = read_profiles($prof_file, 7);

$prof_file = $ENV{'PTPRO_HOURLY'};
my %daily = read_profiles($prof_file, 24);

# open output files
print "Creating output files...\n";
my $output_dir = $ENV{'OUTPUT_DIR'};

# runway files
my $line_loc_fh = open_output("$output_dir/locations/airport_line_locations.csv");
write_line_location_header($line_loc_fh);

my $line_param_fh = open_output("$output_dir/parameters/airport_line_params.csv");
print $line_param_fh "facility_id,facility_name,src_id,src_type,area,fract,relhgt,width,szinit\n";

my $line_tmp_fh = open_output("$output_dir/temporal/airport_line_temporal.csv");
write_temporal_header($line_tmp_fh);

# non-runway files
my $area_loc_fh = open_output("$output_dir/locations/airport_nonrunway_locations.csv");
write_point_location_header($area_loc_fh);

my $area_param_fh = open_output("$output_dir/parameters/airport_nonrunway_params.csv");
print $area_param_fh "facility_id,facility_name,src_id,relhgt,lengthx,lengthy,angle,szinit\n";

my $area_tmp_fh = open_output("$output_dir/temporal/airport_nonrunway_temporal.csv");
write_temporal_header($area_tmp_fh);

# emissions crosswalk file
my $x_fh = open_output("$output_dir/xwalk/airport_srcid_emis.csv");
print $x_fh "state,facility_id,facility_name,smoke_name,ann_value\n";

my %headers;
my @pollutants;
my %records;

print "Processing AERMOD sources...\n";
while (my $line = <$in_fh>) {
  chomp $line;
  next if skip_line($line);
  
  my ($is_header, @data) = parse_report_line($line);

  if ($is_header) {
    parse_header(\@data, \%headers, \@pollutants, 'Plt Name');
    next;
  }
  
  my $state = substr($data[$headers{'Region'}], 7, 2);
  my $plant_id = $data[$headers{'Plant ID'}];
  
  # store all the records by plant ID
  push @{$records{$state}{$plant_id}}, \@data;
}

for my $state (sort keys %records) {
  for my $plant_id (sort keys %{$records{$state}}) {
    my @data = @{$records{$state}{$plant_id}[0]};
  
    # if facility has multiple records, determine which record has most emissions
    if (scalar(@{$records{$state}{$plant_id}}) > 1) {
      my @facility_emissions;  # total emissions for the facility by pollutant
  
      my $selected = 0;
      my $prev_emissions = 0;
      for (my $index = 0; $index < scalar(@{$records{$state}{$plant_id}}); $index++) {
        my $record_emissions = 0;
        foreach my $poll (@pollutants) {
          my $value = $records{$state}{$plant_id}[$index][$headers{$poll}];
          $record_emissions += $value;
          $facility_emissions[$headers{$poll}] += $value;
        }
        if ($record_emissions > $prev_emissions) {
          $selected = $index;
          $prev_emissions = $record_emissions;
        }
      }
    
      @data = @{$records{$state}{$plant_id}[$selected]};
      foreach my $poll (@pollutants) {
        $data[$headers{$poll}] = $facility_emissions[$headers{$poll}];
      }
    }
  
    my $is_runway = exists $runways{$plant_id};
  
    my $src_id;
    if ($is_runway) {
      $src_id = 'AP';
    } else {
      $src_id = 'AP1';
    }
  
    # calculate runway areas
    foreach my $runway (@{$runways{$plant_id}}) {
      my ($zone1, $start_x, $start_y) = latlon_to_utm(23, $runway->{'start_y'}, $runway->{'start_x'});
      my ($zone2, $end_x, $end_y) = latlon_to_utm(23, $runway->{'end_y'}, $runway->{'end_x'});
    
      my $length = sqrt(($end_x - $start_x)**2 + ($end_y - $start_y)**2);
      $runway->{'area'} = $length * $runway->{'width'};
    
      $runway->{'utm_start_x'} = $start_x;
      $runway->{'utm_start_y'} = $start_y;
      $runway->{'utm_end_x'} = $end_x;
      $runway->{'utm_end_y'} = $end_y;
    }

    my @common;
    push @common, $plant_id;
    push @common, '"' . $data[$headers{'Plt Name'}] . '"';
    push @common, $src_id unless $is_runway;

    # prepare location output
    if ($is_runway) {
      my $runway_ct = 1;
      foreach my $runway (@{$runways{$plant_id}}) {
        my @output = $state;
        push @output, @common;
        push @output, $src_id . $runway_ct;
        push @output, $runway->{'utm_start_x'};
        push @output, $runway->{'utm_start_y'};
        push @output, $runway->{'utm_end_x'};
        push @output, $runway->{'utm_end_y'};
        push @output, $data[$headers{'UTM Zone'}];
        push @output, $data[$headers{'X cell'}];
        push @output, $data[$headers{'Y cell'}];
        print $line_loc_fh join(',', @output) . "\n";
        $runway_ct++;
      }
    } else {
      my @output = $state;
      push @output, @common;
      push @output, $data[$headers{'Lambert-X'}];
      push @output, $data[$headers{'Lambert-Y'}];
      push @output, $data[$headers{'Longitude'}];
      push @output, $data[$headers{'Latitude'}];
      push @output, $data[$headers{'UTM_X'}];
      push @output, $data[$headers{'UTM_Y'}];
      push @output, $data[$headers{'UTM Zone'}];
      push @output, $data[$headers{'X cell'}];
      push @output, $data[$headers{'Y cell'}];
      print $area_loc_fh join(',', @output) . "\n";
    }

    # prepare parameters output
    if ($is_runway) {
      my $num_runways = scalar (@{$runways{$plant_id}});
      my $runway_ct = 1;
      foreach my $runway (@{$runways{$plant_id}}) {
        my @output = @common;
        push @output, $src_id . $runway_ct;
        push @output, 'LINE';
        push @output, $runway->{'area'};
        push @output, 1 / $num_runways;
        push @output, 3;
        push @output, $runway->{'width'};
        push @output, 3;
        print $line_param_fh join (',', @output) . "\n";
        $runway_ct++;
      }
    } else {
      my @output = @common;
      push @output, qw(3 10 10 0 3);
      print $area_param_fh join (',', @output) . "\n";
    }

    # prepare temporal profiles output
    my ($qflag, @factors) = get_factors(\%headers, \@data, \%monthly, \%weekly, \%daily);
  
    if ($is_runway) {
      my $runway_ct = 1;
      foreach my $runway (@{$runways{$plant_id}}) {
        my @output = @common;
        push @output, $src_id . $runway_ct;
        push @output, $qflag;
        push @output, map { sprintf('%.4f', $_) } @factors;
        print $line_tmp_fh join(',', @output) . "\n";
      }
    } else {
      my @output = @common;
      push @output, $qflag;
      push @output, map { sprintf('%.4f', $_) } @factors;
      print $area_tmp_fh join(',', @output) . "\n";
    }
  
    # prepare crosswalk output
    pop @common if !$is_runway;
    foreach my $poll (@pollutants) {
      next if $data[$headers{$poll}] == 0.0;
    
      my @output = $state;
      push @output, @common;
      push @output, $poll;
      push @output, $data[$headers{$poll}];
      print $x_fh join(',', @output) . "\n";
    }
  }
}

close $in_fh;
close $line_loc_fh;
close $line_param_fh;
close $line_tmp_fh;
close $area_loc_fh;
close $area_param_fh;
close $area_tmp_fh;
close $x_fh;

print "Done.\n";