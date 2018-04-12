#!/usr/bin/perl

use strict;
use warnings;

use Date::Simple qw(leap_year);
use Geo::Coordinates::UTM qw(latlon_to_utm latlon_to_utm_force_zone);

require 'aermod.subs';
require 'aermod_pt.subs';

my $sector = $ENV{'SECTOR'} || 'ptegu';

my @days_in_month = (0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31);

# check environment variables
foreach my $envvar (qw(REPORT REP_XWALK PHOUR_OUT YEAR 
                       PTPRO_MONTHLY PTPRO_DAILY PTPRO_HOURLY_WINTER PTPRO_HOURLY_SUMMER OUTPUT_DIR)) {
  die "Environment variable '$envvar' must be set" unless $ENV{$envvar};
}

# adjust for leap years
my $year = $ENV{'YEAR'};
if (leap_year($year)) {
  $days_in_month[2] = 29;
}

# open report file
my $in_fh = open_input($ENV{'REPORT'});
my $inx_fh = open_input($ENV{'REP_XWALK'});

# load hourly factors file
print "Reading hourly factors...\n";
my $hour_fh = open_input($ENV{'PHOUR_OUT'});

my %hourly;
while (my $line = <$hour_fh>) {
  chomp $line;
  next unless $line;
  
  my @data = split(/,/, $line);
  my $id = shift @data;
  my $data_year = shift @data;
  die "PHOUR_OUT year ($data_year) doesn't match YEAR setting ($year)" unless $data_year eq $year;
  
  # adjust factors so they sum to 1
  my $sum = sum(@data);
  unless ($sum == 0) {
    @data = map { $_ / $sum } @data;
  }
  
  $id =~ s/^\s+//;
  $id =~ s/\s+$//;
  $hourly{$id} = \@data;
}
close $hour_fh;

# load temporal profiles
print "Reading temporal profiles...\n";
my $prof_file = $ENV{'PTPRO_MONTHLY'};
my %monthly = read_profiles($prof_file, 12);

$prof_file = $ENV{'PTPRO_DAILY'};
my %dayofmonth = read_dom_profiles($prof_file, 31, \@days_in_month);

$prof_file = $ENV{'PTPRO_HOURLY_WINTER'};
my %daily_winter = read_profiles($prof_file, 24);

$prof_file = $ENV{'PTPRO_HOURLY_SUMMER'};
my %daily_summer = read_profiles($prof_file, 24);

# open output files
print "Creating output files...\n";
my $output_dir = $ENV{'OUTPUT_DIR'};

my $loc_fh = open_output("$output_dir/locations/${sector}_location.csv");
write_point_location_header($loc_fh);

my $pt_fh = open_output("$output_dir/parameters/${sector}_point_srcparam.csv");
write_point_srcparam_header($pt_fh);

my $ar_fh = open_output("$output_dir/parameters/${sector}_fug_srcparam.csv");
write_fug_srcparam_header($ar_fh);

my $x_fh = open_output("$output_dir/xwalk/${sector}_srcid_emis.csv");
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

my %headers;
my @pollutants;
my %facility_data; # for each plant_id, a hash with count, utm zone, max emissions, x cell, and y cell
my @locations;
my %hourly_files;

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

  my $src_id = 'SE' . sprintf('%03d', $facility_data{$plant_id}{'count'});
  
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
  unless (exists $hourly_files{$plant_id}) {
    my $fh = open_output("$output_dir/temporal/${plant_id}_${state}_hourly.csv");
    print $fh "facility_id,src_id,year,month,day,hour,hour_factor,stktemp,stkvel\n";
    $hourly_files{$plant_id} = $fh;
  }
  my $fh = $hourly_files{$plant_id};
  
  my $factor_ref;
  my $prof = $data[$headers{'Monthly Prf'}];
  if ($prof =~ /^HR/) {
    $prof =~ s/^HR0*//;
    die "Missing hourly data for source: $prof" unless exists $hourly{$prof};
    $factor_ref = $hourly{$prof};

  } else {
    my $monthly_prof = $data[$headers{'Monthly Prf'}];
    die "Unknown monthly profile code: $monthly_prof" unless exists $monthly{$monthly_prof};
    my @monthly_factors = @{$monthly{$monthly_prof}};
    
    my $dom_prof = $data[$headers{'Day-Month Prf'}];
    die "Unknown day-of-month profile code: $dom_prof" unless exists $dayofmonth{$dom_prof};
    my %dom_factors = %{$dayofmonth{$dom_prof}};
    
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
    die "Invalid hourly profile code: $monday_prof" unless ($monday_prof =~ /w|s$/);
    $monday_prof =~ s/w|s$//;
    die "Unknown hourly profile code: $monday_prof" unless exists $daily_winter{$monday_prof.'w'};
    my @hourly_winter_factors = @{$daily_winter{$monday_prof.'w'}};
    die "Unknown hourly profile code: $monday_prof" unless exists $daily_summer{$monday_prof.'s'};
    my @hourly_summer_factors = @{$daily_summer{$monday_prof.'s'}};

    my @factors;
    my $month = 1;
    foreach my $month_factor (@monthly_factors) {
      # use summer factors for May through Sep
      my @hourly_factors = @hourly_winter_factors;
      @hourly_factors = @hourly_summer_factors if ($month >= 5 && $month <= 9);

      # note: factors average to 1 rather than summing to 1, so adjust when fractions are needed
      my $day_of_month = 1;
      foreach my $dom_factor (@{$dom_factors{$month}}) {
        last if $day_of_month > $days_in_month[$month];
        push @factors, map { $_ * $dom_factor * $month_factor } @hourly_factors;
        $day_of_month++;
      }
      $month++;
    }
    $factor_ref = \@factors;
  }
  
  my $date = Date::Simple->new($year.'-01-01');
  my $hour = 0;
  for my $factor (@$factor_ref) {
    $factor =~ s/^\s+//;
    @output = ($plant_id, $src_id);
    push @output, $date->year;
    push @output, $date->month;
    push @output, $date->day;
    push @output, $hour+1;
    push @output, $factor;
    push @output, $data[$headers{'Stk Tmp'}];
    push @output, $data[$headers{'Stk Vel'}];
    print $fh join (',', @output) . "\n";
    
    $hour++;
    if ($hour == 24) {
      $hour = 0;
      $date = $date->next;
    }
  }
  
  # prepare crosswalk output
  foreach my $poll (@pollutants) {
    @output = $state;
    push @output, @common;
    splice @output, 3, 0, $data[$headers{'Source type'}];
    push @output, $poll;
    push @output, $data[$headers{$poll}];
    print $x_fh join(',', @output) . "\n";
  }

  # prepare inventory source output
  foreach my $src_data (@{$rep_xwalk{$line_num}}) {
    my $xstate = substr(@{$src_data}[2], 7, 2);
    die "Report and crosswalk mismatch at data line $line_num" unless $xstate eq $state && @{$src_data}[3] eq $plant_id;
  
    @output = $state;
    push @output, $plant_id;
    push @output, '"' . $data[$headers{'Plt Name'}] . '"';
    push @output, @{$src_data}[4]; # unit ID
    push @output, @{$src_data}[6]; # process ID
    push @output, @{$src_data}[5]; # release point
    push @output, $src_id;
    print $src_fh join(',', @output) . "\n";
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
close $loc_fh;
close $pt_fh;
close $ar_fh;
close $x_fh;
close $src_fh;

print "Done.\n";
