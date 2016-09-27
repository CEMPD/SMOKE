#!/usr/bin/perl

use strict;
use warnings;

use Scalar::Util qw(looks_like_number);
use List::Util qw(sum);

require 'aermod.subs';

# check environment variables
foreach my $envvar (qw(REPORT PTPRO_MONTHLY PTPRO_WEEKLY PTPRO_HOURLY OUTPUT_DIR)) {
  die "Environment variable '$envvar' must be set" unless $ENV{$envvar};
}

# open report file
my $input = $ENV{'REPORT'};
open (my $in_fh, '<', $input) or die "Could not open file '$input' $!";

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
my $locations = "$output_dir/point_location.csv";
open (my $loc_fh, '>', $locations) or die "Could not open file '$locations' $!";
print $loc_fh "facility_id,facility_name,src_id,grid_x,grid_y,longitude,latitude,utm_x,utm_y,utm_zone,col,row\n";

my $point_params = "$output_dir/point_point_srcparam.csv";
open (my $pt_fh, '>', $point_params) or die "Could not open file '$point_params' $!";
print $pt_fh "facility_id,facility_name,src_id,aermod_src_type,height,temp,velocity,diameter\n";

my $area_params = "$output_dir/point_fug_srcparam.csv";
open (my $ar_fh, '>', $area_params) or die "Could not open file '$area_params' $!";
print $ar_fh "facility_id,facility_name,src_id,aermod_src_type,rel_ht,x_length,y_length,angle,szinit\n";

my $temporal = "$output_dir/point_temporal.csv";
open (my $tmp_fh, '>', $temporal) or die "Could not open file '$temporal' $!";
print $tmp_fh "facility_id,facility_name,src_id,qflag,Scalar1,Scalar2,Scalar3,Scalar4,Scalar5,Scalar6,Scalar7,Scalar8,Scalar9,Scalar10,Scalar11,Scalar12\n";

my $crosswalk = "$output_dir/point_srcid_emis.csv";
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

  my $src_id = 'SN' . $facilities{$plant_id};
  
  my @common;
  push @common, $plant_id;
  push @common, '"' . $data[$headers{'Plt Name'}] . '"';
  push @common, $src_id;
  
  # prepare location output
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
  print $loc_fh join(',', @output) . "\n";
  
  # prepare parameters output
  my $erp_type = $data[$headers{'Emis Release Type'}];
  @output = @common;
  if ($erp_type eq '01') {
    push @output, 'AREA';
    push @output, $data[$headers{'Fug Ht'}];
    push @output, $data[$headers{'Fug Wdt'}];
    push @output, $data[$headers{'Fug Len'}];
    push @output, $data[$headers{'Fug Ang'}];
    
    # calculate szinit
    my $fug_ht = $data[$headers{'Fug Ht'}];
    if ($fug_ht > 10) {
      push @output, sprintf('%.2f', $fug_ht / 4.3);
    } else {
      push @output, 0.00;
    }

    print $ar_fh join(',', @output) . "\n";
  } else {
    if ($erp_type eq '03' ||
        $erp_type eq '06') {
      push @output, 'POINTHOR';
    } elsif ($erp_type eq '05') {
      push @output, 'POINTCAP';
    } else {
      push @output, 'POINT';
    }
    push @output, $data[$headers{'Stk Ht'}];
    push @output, $data[$headers{'Stk Tmp'}];
    push @output, $data[$headers{'Stk Vel'}];
    push @output, $data[$headers{'Stk Dm'}];
    print $pt_fh join(',', @output) . "\n";
  }
  
  # prepare temporal profiles output
  my ($qflag, @factors) = get_factors(\%headers, \@data, \%monthly, \%weekly, \%daily);
  
  @output = @common;
  push @output, $qflag;
  push @output, map { sprintf('%.4f', $_) } @factors;
  print $tmp_fh join(',', @output) . "\n";
  
  # prepare crosswalk output
  foreach my $poll (@pollutants) {
    @output = @common;
    push @output, $poll;
    push @output, $data[$headers{$poll}];
    print $x_fh join(',', @output) . "\n";
  }
}

close $in_fh;
close $loc_fh;
close $pt_fh;
close $ar_fh;
close $tmp_fh;
close $x_fh;

print "Done.\n";
