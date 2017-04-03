#!/usr/bin/perl

use strict;
use warnings;

use Date::Simple ();

require 'aermod.subs';
require 'aermod_pt.subs';

my @days_in_month = (0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31);

# check environment variables
foreach my $envvar (qw(REPORT PHOUR_OUT 
                       PTPRO_MONTHLY PTPRO_DAILY PTPRO_HOURLY_WINTER PTPRO_HOURLY_SUMMER OUTPUT_DIR)) {
  die "Environment variable '$envvar' must be set" unless $ENV{$envvar};
}

# open report file
my $in_fh = open_input($ENV{'REPORT'});

# load hourly factors file
print "Reading hourly factors...\n";
my $hour_fh = open_input($ENV{'PHOUR_OUT'});

my %hourly;
while (my $line = <$hour_fh>) {
  chomp $line;
  next unless $line;
  
  my @data = split(/,/, $line);
  my $id = shift @data;
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

my $loc_fh = open_output("$output_dir/locations/ptegu_location.csv");
write_point_location_header($loc_fh);

my $pt_fh = open_output("$output_dir/parameters/ptegu_point_srcparam.csv");
print $pt_fh "facility_id,facility_name,src_id,aermod_src_type,height,temp,velocity,diameter\n";

my $ar_fh = open_output("$output_dir/parameters/ptegu_fug_srcparam.csv");
print $ar_fh "facility_id,facility_name,src_id,aermod_src_type,rel_ht,x_length,y_length,angle,szinit\n";

my $x_fh = open_output("$output_dir/xwalk/ptegu_srcid_emis.csv");
write_crosswalk_header($x_fh);

my %headers;
my @pollutants;
my %facilities;
my %hourly_files;

print "Processing AERMOD sources...\n";
while (my $line = <$in_fh>) {
  chomp $line;
  next if skip_line($line);

  my ($is_header, @data) = parse_report_line($line);
  
  if ($is_header) {
    parse_header(\@data, \%headers, \@pollutants, 'Plt Name');
    next;
  }
  
  # add to source count for current facility
  my $plant_id = $data[$headers{'Plant ID'}];
  unless (exists $facilities{$plant_id}) {
    $facilities{$plant_id} = 0;
  }
  $facilities{$plant_id}++;

  my $src_id = 'SE' . $facilities{$plant_id};
  
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
  unless (exists $hourly_files{$plant_id}) {
    my $region = $data[$headers{'Region'}];
    $region = substr($region, -6);
    my $state = substr($region, 1, 2);
    my $fh = open_output("$output_dir/temporal/${plant_id}_${state}_hourly.csv");
    print $fh "facility_id,src_id,year,month,day,hour,hour_factor,stktemp,stkvel\n";
    $hourly_files{$plant_id} = $fh;
  }
  my $fh = $hourly_files{$plant_id};
  
  my $factor_ref;
  my $prof = $data[$headers{'Monthly Prf'}];
  my $year = 2014;
  if ($prof =~ /^HR/) {
    $prof =~ s/^HR0*//;
    die "Missing hourly data for source: $prof" unless exists $hourly{$prof};
    $factor_ref = $hourly{$prof};
    $year = shift @$factor_ref;

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
    push @output, $hour;
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
close $x_fh;

print "Done.\n";
