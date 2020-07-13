#!/usr/bin/perl

use strict;
use warnings;

use Text::CSV ();
use Geo::Coordinates::UTM qw(latlon_to_utm_force_zone);

require 'aermod.subs';
require 'aermod_pt.subs';

# check environment variables
foreach my $envvar (qw(OUTPUT_DIR)) {
  die "Environment variable '$envvar' must be set" unless $ENV{$envvar};
}

my $output_dir = $ENV{'OUTPUT_DIR'};

# open sector-specific files
my $ptegu_loc_fh = open_input("$output_dir/locations/ptegu_location.csv");
my $ptnon_loc_fh = open_input("$output_dir/locations/ptnonipm_location.csv");

my $ptegu_pt_fh = open_input("$output_dir/parameters/ptegu_point_srcparam.csv");
my $ptnon_pt_fh = open_input("$output_dir/parameters/ptnonipm_point_srcparam.csv");

my $ptegu_ar_fh = open_input("$output_dir/parameters/ptegu_fug_srcparam.csv");
my $ptnon_ar_fh = open_input("$output_dir/parameters/ptnonipm_fug_srcparam.csv");

my $ptegu_x_fh = open_input("$output_dir/xwalk/ptegu_process_releasept_emis.csv");
my $ptnon_x_fh = open_input("$output_dir/xwalk/ptnonipm_process_releasept_emis.csv");

# open combined output files
print "Creating output files...\n";
my $loc_fh = open_output("$output_dir/locations/point_location.csv");
write_point_location_header($loc_fh);

my $pt_fh = open_output("$output_dir/parameters/point_point_srcparam.csv");
write_point_srcparam_header($pt_fh);

my $ar_fh = open_output("$output_dir/parameters/point_fug_srcparam.csv");
write_fug_srcparam_header($ar_fh);

my $x_fh = open_output("$output_dir/xwalk/point_combined_process_releasept_emis.csv");
write_crosswalk_header($x_fh);

# read ptegu data
my $csv_parser = Text::CSV->new();
$csv_parser->eol("\n");

my $header = $csv_parser->getline($ptegu_loc_fh);
$csv_parser->column_names(@$header);

my %ptegu;
while (my $row = $csv_parser->getline_hr($ptegu_loc_fh)) {
  my $plant_id = $row->{'facility_id'};
  unless (exists $ptegu{$plant_id}) {
    $ptegu{$plant_id} = $row;
  }
  
  # output line to combined file (can't use Text::CSV->print_hr due to requirement
  # that facility name always be in double quotes)
  my @output;
  foreach my $name ($csv_parser->column_names) {
    if ($name eq 'facility_name') {
      push @output, qq("$row->{$name}");
    } else {
      push @output, $row->{$name};
    }
  }
  print $loc_fh join(',', @output) . "\n";
}

# read ptnonipm data
$csv_parser->getline($ptnon_loc_fh); # skip header line
while (my $row = $csv_parser->getline_hr($ptnon_loc_fh)) {
  my $plant_id = $row->{'facility_id'};
  if (exists $ptegu{$plant_id}) {
    # check UTM zone and recalculate utm_x and utm_y if needed
    if ($row->{'utm_zone'} != $ptegu{$plant_id}->{'utm_zone'}) {
      my ($zone, $utm_x, $utm_y) = latlon_to_utm_force_zone(23, $ptegu{$plant_id}->{'utm_zone'}, $row->{'latitude'}, $row->{'longitude'});
      $row->{'utm_x'} = sprintf('%.2f', $utm_x);
      $row->{'utm_y'} = sprintf('%.2f', $utm_y);
      $row->{'utm_zone'} = $ptegu{$plant_id}->{'utm_zone'};
    }
    
    # set grid cell based on ptegu data
    $row->{'col'} = $ptegu{$plant_id}->{'col'};
    $row->{'row'} = $ptegu{$plant_id}->{'row'};
  }
  
  my @output;
  foreach my $name ($csv_parser->column_names) {
    if ($name eq 'facility_name') {
      push @output, qq("$row->{$name}");
    } else {
      push @output, $row->{$name};
    }
  }
  print $loc_fh join(',', @output) . "\n";
}

# combine other files
while (my $line = <$ptegu_pt_fh>) {
  next if $line =~ /^facility_id/; # skip header line
  print $pt_fh $line;
}

while (my $line = <$ptnon_pt_fh>) {
  next if $line =~ /^facility_id/; # skip header line
  print $pt_fh $line;
}

while (my $line = <$ptegu_ar_fh>) {
  next if $line =~ /^facility_id/; # skip header line
  print $ar_fh $line;
}

while (my $line = <$ptnon_ar_fh>) {
  next if $line =~ /^facility_id/; # skip header line
  print $ar_fh $line;
}

while (my $line = <$ptegu_x_fh>) {
  next if $line =~ /^state/; # skip header line
  print $x_fh $line;
}

while (my $line = <$ptnon_x_fh>) {
  next if $line =~ /^state/; # skip header line
  print $x_fh $line;
}

close $ptegu_loc_fh;
close $ptnon_loc_fh;
close $ptegu_pt_fh;
close $ptnon_pt_fh;
close $ptegu_ar_fh;
close $ptnon_ar_fh;
close $ptegu_x_fh;
close $ptnon_x_fh;
close $loc_fh;
close $pt_fh;
close $ar_fh;
close $x_fh;

print "Done.\n";
