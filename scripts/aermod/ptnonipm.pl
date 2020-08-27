#!/usr/bin/perl

use strict;
use warnings;

use Geo::Coordinates::UTM qw(latlon_to_utm latlon_to_utm_force_zone);

require 'aermod.subs';
require 'aermod_pt.subs';

my $sector = $ENV{'SECTOR'} || 'ptnonipm';

# check environment variables
foreach my $envvar (qw(REPORT REP_XWALK REP_SRC PTPRO_MONTHLY PTPRO_WEEKLY PTPRO_HOURLY OUTPUT_DIR)) {
  die "Environment variable '$envvar' must be set" unless $ENV{$envvar};
}

# open report files
my $in_fh = open_input($ENV{'REPORT'});
my $inx_fh = open_input($ENV{'REP_XWALK'});
my $ins_fh = open_input($ENV{'REP_SRC'});

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

my $loc_fh = open_output("$output_dir/locations/${sector}_location.csv");
write_point_location_header($loc_fh);

my $pt_fh = open_output("$output_dir/parameters/${sector}_point_srcparam.csv");
write_point_srcparam_header($pt_fh);

my $ar_fh = open_output("$output_dir/parameters/${sector}_fug_srcparam.csv");
write_fug_srcparam_header($ar_fh);

my $tmp_fh = open_output("$output_dir/temporal/${sector}_temporal.csv");
write_temporal_header($tmp_fh);

my $x_fh = open_output("$output_dir/xwalk/${sector}_process_releasept_emis.csv");
write_crosswalk_header($x_fh);

my $src_fh = open_output("$output_dir/xwalk/${sector}_srcid_xwalk.csv");
write_source_header($src_fh);

my %rep_xwalk;

print "Reading AERMOD-to-inventory source info...\n";
while (my $line = <$inx_fh>) {
  chomp $line;
  $line =~ s/^\s+//;
  $line =~ s/\s+$//;

  next if $line =~ /^GROUPID/;

  my @data = split(/\s+/, $line);
  push @{$rep_xwalk{$data[0]}}, \@data;
}

my %src_headers;
my @src_pollutants;
my %rep_src;
print "Reading inventory source report...\n";
while (my $line = <$ins_fh>) {
  chomp $line;
  next if skip_line($line);
  
  my ($is_header, @data) = parse_report_line($line);
  
  if ($is_header) {
    parse_header(\@data, \%src_headers, \@src_pollutants, 'Char 3');
    next;
  }
  
  my $smoke_id = $data[$src_headers{'Source ID'}];
  $rep_src{$smoke_id} = \@data;
}

my %headers;
my @pollutants;
my %facility_data; # for each plant_id, a hash with count, utm zone, max emissions, x cell, and y cell
my @locations;

print "Processing AERMOD sources...\n";
my $line_num = 0;
while (my $line = <$in_fh>) {
  chomp $line;
  next if skip_line($line);

  my ($is_header, @data) = parse_report_line($line);
  
  if ($is_header) {
    parse_header(\@data, \%headers, \@pollutants, 'Plt Name');
    next;
  }
  
  $line_num++;
  
  # sum emissions and check if all emissions are zero
  my $record_emissions = 0;
  foreach my $poll (@pollutants) {
    $record_emissions += $data[$headers{$poll}];
  }
  next if $record_emissions == 0;
  
  # initialize facility data if needed
  my $plant_id = $data[$headers{'Plant ID'}];
  unless (exists $facility_data{$plant_id}) {
    $facility_data{$plant_id}{'count'} = 0;
    $facility_data{$plant_id}{'max_emissions'} = $record_emissions;
    $facility_data{$plant_id}{'x cell'} = $data[$headers{'X cell'}];
    $facility_data{$plant_id}{'y cell'} = $data[$headers{'Y cell'}];
    my ($zone, $utm_x, $utm_y) = latlon_to_utm(23, $data[$headers{'Latitude'}], $data[$headers{'Longitude'}]);
    $facility_data{$plant_id}{'utm zone'} = $zone;
  }

  # add to source count for current facility
  $facility_data{$plant_id}{'count'}++;
  
  # update grid cell if current record has more emissions
  if ($record_emissions > $facility_data{$plant_id}{'max_emissions'}) {
    $facility_data{$plant_id}{'max_emissions'} = $record_emissions;
    $facility_data{$plant_id}{'x cell'} = $data[$headers{'X cell'}];
    $facility_data{$plant_id}{'y cell'} = $data[$headers{'Y cell'}];
  }

  my $src_id = 'SN' . sprintf('%03d', $facility_data{$plant_id}{'count'});
  
  my @common;
  push @common, $plant_id;
  push @common, '"' . $data[$headers{'Plt Name'}] . '"';
  push @common, $src_id;
  
  # prepare location output
  my $state = substr($data[$headers{'Region'}], 7, 2);
  my @output = $state;
  push @output, @common;
  push @output, $data[$headers{'Lambert-X'}];
  push @output, $data[$headers{'Lambert-Y'}];
  push @output, $data[$headers{'Longitude'}];
  push @output, $data[$headers{'Latitude'}];
  my ($zone, $utm_x, $utm_y) = latlon_to_utm_force_zone(23, $facility_data{$plant_id}{'utm zone'}, $data[$headers{'Latitude'}], $data[$headers{'Longitude'}]);
  my $outzone = $zone;
  $outzone =~ s/\D//g; # strip latitude band designation from UTM zone
  push @output, sprintf('%.2f', $utm_x);
  push @output, sprintf('%.2f', $utm_y);
  push @output, $outzone;
  # grid cell row and column will be added just before output
  push @locations, [ @output ];
  
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
      push @output, sprintf('%.5f', $fug_ht / 4.3);
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
  push @output, map { sprintf('%.8f', $_) } @factors;
  print $tmp_fh join(',', @output) . "\n";
  
  # prepare crosswalk output
  foreach my $xwalk_arrayref (@{$rep_xwalk{$line_num}}) {
    my @xwalk_data = @{$xwalk_arrayref};
    my $xstate = substr($xwalk_data[2], 7, 2);
    die "Report and crosswalk mismatch at data line $line_num" unless $xstate eq $state && $xwalk_data[3] eq $plant_id;
    
    my $smoke_id = $xwalk_data[1];
    die "Missing data for Source ID: $smoke_id in REP_SRC" unless exists $rep_src{$smoke_id};
    my @src_data = @{$rep_src{$smoke_id}};
    
    # output inventory source file for QA
    @output = $state;
    push @output, $plant_id;
    push @output, '"' . $data[$headers{'Plt Name'}] . '"';
    push @output, $xwalk_data[4]; # unit ID
    push @output, $xwalk_data[6]; # process ID
    push @output, $xwalk_data[5]; # release point
    push @output, $src_id;
    print $src_fh join(',', @output) . "\n";
    
    foreach my $poll (@src_pollutants) {
      next if $src_data[$src_headers{$poll}] == 0.0;
      @output = $state;
      push @output, $plant_id;
      push @output, '"' . $data[$headers{'Plt Name'}] . '"';
      push @output, $data[$headers{'Source type'}];
      push @output, $src_id;
      push @output, $xwalk_data[4]; # unit ID
      push @output, $xwalk_data[6]; # process ID
      push @output, $xwalk_data[5]; # release point
      push @output, $poll;
      push @output, $src_data[$src_headers{$poll}];
      print $x_fh join(',', @output) . "\n";
    }
  }
}

# output location records
foreach my $out_ref (@locations) {
  my $plant_id = $out_ref->[1];
  push @$out_ref, $facility_data{$plant_id}{'x cell'};
  push @$out_ref, $facility_data{$plant_id}{'y cell'};
  print $loc_fh join(',', @$out_ref) . "\n";
}

close $in_fh;
close $inx_fh;
close $ins_fh;
close $loc_fh;
close $pt_fh;
close $ar_fh;
close $tmp_fh;
close $x_fh;
close $src_fh;

print "Done.\n";
