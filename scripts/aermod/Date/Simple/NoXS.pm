# Date::Simple::NoXS - used internally by Date::Simple.

use strict;

package Date::Simple;

sub _ymd {
    my ( $class, @args ) = @_;
    my $days = ymd_to_days(@args);
    return unless defined($days);
    return ( bless \$days, $class );
}

sub _d8 {
    my ( $o, $d8 ) = @_;
    my @ymd = $d8 =~ m/^(\d{4})(\d\d)(\d\d)$/ or return undef;
    return $o->_ymd(@ymd);
}

# Precise integer arithmetic functions unfortunately missing from
# Perl's core:

sub _divmod {
    my ( $quot, $int );

    $quot = $_[0] / $_[1];
    $int  = int($quot);
    $int -= 1 if $int > $quot;
    $_[0] %= $_[1];
    return $int;
}

sub _div {
    my ( $quot, $int );

    $quot = $_[0] / $_[1];
    $int  = int($quot);
    return $int - 1 if $int > $quot;
    return $int;
}

sub leap_year {
    my $y = shift;
    return ( ( $y % 4 == 0 ) and ( $y % 400 == 0 or $y % 100 != 0 ) ) || 0;
}

my @days_in_month = (
    [ 0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ],
    [ 0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ],
);

sub days_in_month ($$) {
    my ( $y, $m ) = @_;
    return $days_in_month[ leap_year($y) ][$m];
}

sub validate ($$$) {
    my ( $y, $m, $d ) = @_;

    # any +ve integral year is valid
    return 0 if $y != abs int $y;
    return 0 unless 1 <= $m and $m <= 12;
    return 0 unless 1 <= $d and $d <= $days_in_month[ leap_year($y) ][$m];
    return 1;
}

# Given a year, month, and day, return the canonical day number.
# That is the number of days since 1 January 1970, negative if earlier.
sub ymd_to_days {
    my ( $Y, $M, $D ) = @_;
    my ( $days, $x );

    if (   $M < 1
        || $M > 12
        || $D < 1
        || ( $D > 28 && $D > days_in_month( $Y, $M ) ) ) {
        return undef;
    }

    $days = $D +
      ( undef, -1, 30, 58, 89, 119, 150, 180, 211, 242, 272, 303, 333 )[$M];
    $days += 365 * ( $Y - 1970 );
    $x = ( $M <= 2 ? $Y - 1 : $Y );
    $days += _div( ( $x - 1968 ), 4 );
    $days -= _div( ( $x - 1900 ), 100 );
    $days += _div( ( $x - 1600 ), 400 );
    return $days;
}

sub days_since_1970 { ${ $_[0] } }

# Given a canonical day number (days since 1 Jan 1970), return the
# year, month, and day.
sub days_to_ymd {
    my ($days) = @_;
    my ( $year, $mnum, $mday, $tmp );

    # Shift frame of reference from 1 Jan 1970 to (the imaginary) 1 Mar 0AD.
    $tmp = $days + 719468;

    # Do the math.
    $year = 400 * _divmod( $tmp, 146097 );
    if ( $tmp == 146096 ) {

        # Handle 29 Feb 2000, 2400, ...
        $year += 400;
        $mnum = 2;
        $mday = 29;
    }
    else {
        $year += 100 * _divmod( $tmp, 36524 );
        $year += 4 * _divmod( $tmp, 1461 );
        if ( $tmp == 1460 ) {
            $year += 4;
            $mnum = 2;
            $mday = 29;
        }
        else {
            $year += _divmod( $tmp, 365 );
            $mnum = _divmod( $tmp, 31 );
            $mday = $tmp + ( 1, 1, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5 )[$mnum];
            $tmp = ( 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 31, 28 )[$mnum];
            if ( $mday > $tmp ) {
                $mday -= $tmp;
                $mnum += 1;
            }
            if ( $mnum > 9 ) {
                $mnum -= 9;
                $year += 1;
            }
            else {
                $mnum += 3;
            }
        }
    }
    return ( $year, $mnum, $mday );
}

sub as_ymd { return days_to_ymd( ${ $_[0] } ); }
sub as_d8  { return sprintf( "%04d%02d%02d", &as_ymd ); }
sub as_iso { return sprintf( "%04d-%02d-%02d", &as_ymd ); }

sub year  { return (&as_ymd)[0]; }
sub month { return (&as_ymd)[1]; }
sub day   { return (&as_ymd)[2]; }

sub day_of_week {
    return ( ( ${ $_[0] } + 4 ) % 7 );
}

#------------------------------------------------------------------------------
# the following methods are called by the overloaded operators, so they should
# not normally be called directly.
#------------------------------------------------------------------------------

sub _add {
    my ( $date, $diff ) = @_;

    if ( $diff !~ /^-?\d+$/ ) {
        Carp::croak("Date interval must be an integer");
    }
    my $new_date = bless( \( $$date + $diff ), ref($date) );
    $new_date->default_format( $date->default_format );
    return $new_date;
}

sub _subtract {
    my ( $left, $right, $reverse ) = @_;

    my $new_date;

    if ($reverse) {
        Carp::croak("Can't subtract a date from a non-date");
    }
    if ( ref($right) eq '' && $right =~ /^-?\d+$/ ) {
        $new_date = bless( \( $$left - $right ), ref($left) );
        $new_date->default_format( $left->default_format );
        return $new_date;
    }
    return ( $$left - $$right );
}

sub _compare {
    my ( $left, $right, $reverse ) = @_;

    $right = $left->new($right) || _inval( $left, $right );
    return ( $reverse ? $$right <=> $$left : $$left <=> $$right );
}

sub _eq {
    my ( $left, $right ) = @_;
    return ( ( $right = $left->_new($right) ) && $$right == $$left );
}

sub _ne {
    return ( !&_eq );
}

1;

=head1 NAME

Date::Simple::NoXS - Pure Perl support for Date::Simple.

=head1 SYNOPSIS

    use Date::Simple;

=head1 DESCRIPTION

Used internally by Date::Simple.

=head1 SEE ALSO

L<Date::Simple>.

=cut
