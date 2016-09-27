#!/usr/bin/perl

use strict;
use warnings;

use Scalar::Util qw(looks_like_number);
use List::Util qw(sum);
use Text::CSV;
use Geo::Coordinates::UTM;

require 'aermod.subs';

# check environment variables
foreach my $envvar (qw(REPORT RUNWAYS PTPRO_MONTHLY PTPRO_WEEKLY PTPRO_HOURLY OUTPUT_DIR)) {
  die "Environment variable '$envvar' must be set" unless $ENV{$envvar};
}

# open report file
my $input = $ENV{'REPORT'};
open (my $in_fh, '<', $input) or die "Could not open file '$input' $!";

# load runway data
print "Reading runway data...\n";
my $runways_file = $ENV{'RUNWAYS'};
open (my $runway_fh, '<', $runways_file) or die "Could not open file '$runways_file' $!";

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
my $line_locations = "$output_dir/airport_line_locations.csv";
open (my $line_loc_fh, '>', $line_locations) or die "Could not open file '$line_locations' $!";
print $line_loc_fh "facility_id,facility_name,src_id,xs1,ys1,xs2,ys2,utm_zone,col,row\n";

my $line_params = "$output_dir/airport_line_params.csv";
open (my $line_param_fh, '>', $line_params) or die "Could not open file '$line_params' $!";
print $line_param_fh "facility_id,facility_name,src_id,src_type,area,fract,relhgt,width,szinit\n";

my $line_temporal = "$output_dir/airport_line_temporal.csv";
open (my $line_tmp_fh, '>', $line_temporal) or die "Could not open file '$line_temporal' $!";
print $line_tmp_fh "facility_id,facility_name,src_id,qflag,Scalar1,Scalar2,Scalar3,Scalar4,Scalar5,Scalar6,Scalar7,
Scalar8,Scalar9,Scalar10,Scalar11,Scalar12\n";

# non-runway files
my $area_locations = "$output_dir/airport_nonrunway_locations.csv";
open (my $area_loc_fh, '>', $area_locations) or die "Could not open file '$area_locations' $!";
print $area_loc_fh "facility_id,facility_name,src_id,grid_x,grid_y,longitude,latitude,utm_x,utm_y,utm_zone,col,row\n";

my $area_params = "$output_dir/airport_nonrunway_params.csv";
open (my $area_param_fh, '>', $area_params) or die "Could not open file '$area_params' $!";
print $area_param_fh "facility_id,facility_name,src_id,relhgt,lengthx,lengthy,angle,szinit\n";

my $area_temporal = "$output_dir/airport_nonrunway_temporal.csv";
open (my $area_tmp_fh, '>', $area_temporal) or die "Could not open file '$area_temporal' $!";
print $area_tmp_fh "facility_id,facility_name,src_id,qflag,Scalar1,Scalar2,Scalar3,Scalar4,Scalar5,Scalar6,Scalar7,
Scalar8,Scalar9,Scalar10,Scalar11,Scalar12\n";

# emissions crosswalk file
my $crosswalk = "$output_dir/airport_srcid_emis.csv";
open (my $x_fh, '>', $crosswalk) or die "Could not open file '$crosswalk' $!";
print $x_fh "facility_id,rel_point_id,scc,src_id,poll,ann_value\n";

my %headers;
my @pollutants;
my %facilities;

print "Processing AERMOD sources...\n";
while (my $line = <$in_fh>) {
  chomp $line;
  next if skip_line($line);
  
  my ($is_header, @data) = parse_report_line($line);

  if ($is_header) {
    parse_header(\@data, \%headers, \@pollutants);
    next;
  }

  # add to source count for current facility
  my $plant_id = $data[$headers{'Plant ID'}];
  unless (exists $facilities{$plant_id}) {
    $facilities{$plant_id} = 0;
  }
  $facilities{$plant_id}++;
  
  my $is_runway = exists $runways{$plant_id};
  
  my $src_id;
  if ($is_runway) {
    $src_id = 'S';
  } else {
    $src_id = 'AP' . $facilities{$plant_id};
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
    my $runway_ct = 0;
    foreach my $runway (@{$runways{$plant_id}}) {
      my @output = @common;
      push @output, $src_id . ($facilities{$plant_id} + $runway_ct);
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
    my @output = @common;
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
    my $runway_ct = 0;
    foreach my $runway (@{$runways{$plant_id}}) {
      my @output = @common;
      push @output, $src_id . ($facilities{$plant_id} + $runway_ct);
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
    push @output, qw(3 10 10 0);
    print $area_param_fh join (',', @output) . "\n";
  }

  # prepare temporal profiles output
  my ($qflag, @factors) = get_factors(\%headers, \@data, \%monthly, \%weekly, \%daily);
  
  if ($is_runway) {
    my $runway_ct = 0;
    foreach my $runway (@{$runways{$plant_id}}) {
      my @output = @common;
      push @output, $src_id . ($facilities{$plant_id} + $runway_ct);
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
  if ($is_runway) {
    my $runway_ct = 0;
    foreach my $runway (@{$runways{$plant_id}}) {
      foreach my $poll (@pollutants) {
        my @output = @common;
        push @output, $src_id . ($facilities{$plant_id} + $runway_ct);
        push @output, $poll;
        push @output, $data[$headers{$poll}];
        print $x_fh join(',', @output) . "\n";
      }
    }
  } else {
    foreach my $poll (@pollutants) {
      my @output = @common;
      push @output, $poll;
      push @output, $data[$headers{$poll}];
      print $x_fh join(',', @output) . "\n";
    }
  }
  
  # update count of sources at facility for individual runways
  if ($is_runway) {
    $facilities{$plant_id} += scalar (@{$runways{$plant_id}});
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
