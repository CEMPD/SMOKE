# Date::Simple - a simple date object

package Date::Simple;

BEGIN {
    $VERSION = '3.03';
}

use Exporter ();
@ISA = ('Exporter');

@EXPORT_OK = qw(
  today ymd d8 leap_year days_in_month
  date date_fmt date_d8 date_iso
);

%EXPORT_TAGS = ( 'all' => \@EXPORT_OK );

# Try to load the C code.  If that fails, fall back to Date::Simple::NoXS.
if ( !defined(&_add) ) {
    my $err = $Date::Simple::NoXS;
    unless ($err) {

        # Use DynaLoader instead of XSLoader for pre-5.005.
        local ($@);
        local @ISA = ('DynaLoader');
        require DynaLoader;
        eval { __PACKAGE__->bootstrap($VERSION); };
        $err = $@;
    }
    if ($err) {
        $Date::Simple::NoXs = 1;
        require Date::Simple::NoXS;
    }
}

use strict;
use Carp ();
use overload
  '+'    => '_add',
  '-'    => '_subtract',
  '=='   => '_eq',
  '!='   => '_ne',
  '<=>'  => '_compare',
  'eq'   => '_eq',
  'ne'   => '_ne',
  'cmp'  => '_compare',
  'bool' => sub { 1 },
  '""'   => 'as_iso';

use Scalar::Util qw(refaddr reftype);
use warnings::register;
require Date::Simple::Fmt;
require Date::Simple::ISO;
require Date::Simple::D8;

sub d8 {

    # called as function
    if ( $#_ == 0 ) {
        return __PACKAGE__->_d8(@_);
    }

    # called as method
    else {
        if ( ref $_[0] eq 'SCALAR' ) {
            return $_[0]->SUPER::_d8(@_);
        }
        else {
            return $_[0]->_d8(@_);
        }
    }
}

sub today {
    if ( $#_ == -1 ) {
        return __PACKAGE__->_today(@_);
    }
    else {
        return shift->_today(@_);
    }
}

sub ymd {

    # called as function
    if ( $#_ == 2 ) {
        return __PACKAGE__->_ymd(@_);
    }

    # called as method
    else {
        if ( ref $_[0] eq 'SCALAR' ) {
            return $_[0]->SUPER::_ymd(@_);
        }
        else {
            return $_[0]->_ymd(@_);
        }
    }
}

sub _today {
    my ( $y, $m, $d ) = (localtime)[ 5, 4, 3 ];
    $y += 1900;
    $m += 1;
    return $_[0]->_ymd( $y, $m, $d );
}

sub _inval {
    my ($first);
    $first = shift;
    Carp::croak( "Invalid "
          . ( ref($first) || $first )
          . " constructor args: ('"
          . join( "', '", @_ )
          . "')" );
}

sub _new {
    my ( $that, @ymd ) = @_;

    my $class = ref($that) || $that;

    if ( @ymd == 1 ) {
        my $x = $ymd[0];
        if ( ref $x and reftype($x) eq 'ARRAY' ) {
            @ymd = @$x;
        }
        elsif ( UNIVERSAL::isa( $x, __PACKAGE__ ) ) {
            return ($x);
        }
        elsif ($x =~ /^(\d\d\d\d)-(\d\d)-(\d\d)$/
            || $x =~ /^(\d\d\d\d)(\d\d)(\d\d)$/ ) {
            @ymd = ( $1, $2, $3 );
        }
        else {
            return (undef);
        }
    }    # we fall through here...

    # note we can end up here is they pass in [] as the date
    return $class->_today() unless @ymd;

    # to get here, we had one arg which was split,
    # or 3 in the first place
    if ( @ymd == 3 ) {
        my $days = ymd_to_days(@ymd);
        return undef if !defined($days);
        return ( bless( \$days, $class ) );
    }

    $class->_inval(@ymd);
}

sub date { scalar __PACKAGE__->_new(@_) }

sub date_fmt {
    my $format = shift;
    my $obj    = Date::Simple::Fmt->_new(@_);
    $obj->default_format($format)
      if $obj;
    $obj;
}

sub date_d8  { scalar Date::Simple::D8->_new(@_) }
sub date_iso { scalar Date::Simple::ISO->_new(@_) }

# Same as date() but it's a method and croaks on error if called with
# one arg.
sub new {
    my ( $class, $date );

    $date = &_new;
    if ( !$date && scalar(@_) == 1 ) {
        Carp::croak( "'" . shift() . "' is not a valid ISO formated date" );
    }
    return ($date);
}

sub next { return ( $_[0] + 1 ); }
sub prev { return ( $_[0] - 1 ); }

sub _gmtime {
    my ( $y, $m, $d ) = days_to_ymd( ${ $_[0] } );
    $y -= 1900;
    $m -= 1;
    return ( 0, 0, 0, $d, $m, $y );
}

BEGIN {
    our $Standard_Format = "%Y-%m-%d";
    my %fmts = (    # Inside out parameter
        'Date::Simple'      => $Standard_Format,
        'Date::Simple::ISO' => $Standard_Format,
        'Date::Simple::D8'  => "%Y%m%d",
        'Date::Simple::Fmt' => $Standard_Format,
    );

    sub format {
        my ( $self, $format ) = @_;

        $format =
             $fmts{ refaddr($self) || '' }
          || $fmts{ ref($self) }
          || $Standard_Format
          if @_ == 1;

        return "$self" unless defined($format);
        require POSIX;
        local $ENV{TZ} = 'UTC+0';
        return POSIX::strftime( $format, _gmtime($self) );
    }

    sub strftime { &format }
    sub as_str   { &format }

    sub default_format {
        my ( $self, $val ) = @_;

        my $o = refaddr($self) || $self;

        if ( @_ > 1 ) {
            $fmts{$o} = $val;
            warnings::warnif "Setting class specific date format '$o' to" . "'"
              . ( defined $val ? $val : 'undef' ) . "'"
              unless ref $self;
        }

        return $fmts{$o} || $Standard_Format;
    }

    sub DESTROY {
        delete $fmts{ refaddr $_[0] };
    }
}

1;

=head1 NAME

Date::Simple - a simple date object

=head1 SYNOPSIS

    use Date::Simple ('date', 'today');

    # Difference in days between two dates:
    $diff = date('2001-08-27') - date('1977-10-05');

    # Offset $n days from now:
    $date = today() + $n;
    print "$date\n";  # uses ISO 8601 format (YYYY-MM-DD)

    use Date::Simple ();
    my $date  = Date::Simple->new('1972-01-17');
    my $year  = $date->year;
    my $month = $date->month;
    my $day   = $date->day;

    use Date::Simple (':all');
    my $date2 = ymd($year, $month, $day);
    my $date3 = d8('19871218');
    my $today = today();
    my $tomorrow = $today + 1;
    if ($tomorrow->year != $today->year) {
        print "Today is New Year's Eve!\n";
    }

    if ($today > $tomorrow) {
        die "warp in space-time continuum";
    }

    print "Today is ";
    print(('Sun','Mon','Tues','Wednes','Thurs','Fri','Satur')
          [$today->day_of_week]);
    print "day.\n";

    # you can also do this:
    ($date cmp "2001-07-01")
    # and this
    ($date <=> [2001, 7, 1])

=begin text

INSTALLATION

 If your system has the "make" program or a clone:

     perl Makefile.PL
     make
     make test
     make install

 If you lack "make", copy the "lib/Date" directory to your module
 directory (run "perl -V:sitelib" to find it).

 If "make test" fails, perhaps it means your system can't compile C
 code.  Try:

     make distclean
     perl Makefile.PL noxs
     make
     make test
     make install

 This will use the pure-Perl implementation.

=end text

=head1 DESCRIPTION

Dates are complex enough without times and timezones.  This module may
be used to create simple date objects.  It handles:

=over 4

=item Validation.

Reject 1999-02-29 but accept 2000-02-29.

=item Interval arithmetic.

How many days were between two given dates?  What date comes N days
after today?

=item Day-of-week calculation.

What day of the week is a given date?

=item Transparent date formatting.

How should a date object be formatted.

=back

It does B<not> deal with hours, minutes, seconds, and time zones.

A date is uniquely identified by year, month, and day integers within
valid ranges.  This module will not allow the creation of objects for
invalid dates.  Attempting to create an invalid date will return
undef.  Month numbering starts at 1 for January, unlike in C and Java.
Years are 4-digit.

Gregorian dates up to year 9999 are handled correctly, but we rely on
Perl's builtin C<localtime> function when the current date is
requested.  On some platforms, C<localtime> may be vulnerable to
rollovers such as the Unix C<time_t> wraparound of 18 January 2038.

Overloading is used so you can compare or subtract two dates using
standard numeric operators such as C<==>, and the sum of a date object
and an integer is another date object.

Date::Simple objects are immutable.  After assigning C<$date1> to
C<$date2>, no change to C<$date1> can affect C<$date2>.  This means,
for example, that there is nothing like a C<set_year> operation, and
C<$date++> assigns a new object to C<$date>.

This module contains various undocumented functions.  They may not be
available on all platforms and are likely to change or disappear in
future releases.  Please let the author know if you think any of them
should be public.

=head2 Controlling output format.

As of version 3.0 new ways of controlling the output formats of Date::Simple
objects has been provided. However Date::Simple has traditionally provided
few ways of stringification, a primary one via the format() method and another
primary one via direct stringification. However the later is currently
implemented as an XS routine and the former is implemented through a perl routine.
This means that using format() is more expensive than stringification and
that the stringification format is class specific.

In order to alleviate some of these problems a new mechanism has been introduced
to Date::Simple that allows for a per object level format default. In addition
a set of utility classes that have different stringification overloads provided.
These classes are simple subclasses of Date::Simple and beside the default format()
and the overloaded stringification behaviour are identical to Date::Simple. In fact
one is totally identical to Date::Simple and is provided mostly for completeness.

The classes included are:

=over 4

=item Date::Simple::ISO

Identical to Date::Simple in every respect but name.

=item Date::Simple::D8

Uses the D8 format (%Y%m%d) as the default format for printing. Uses XS for the
overloaded stringification.

=item Date::Simple::Fmt

Uses the perl implemented format() as the default stringification mechanism. The first
argument to the constructor is expected to be the format to use for the object.

=back

B<NOTE> its important to remember that the primary difference between the behaviour
of objects of the different classes is how they are stringified when quoted, and what
date format is used by default when the format() method is called. Nothing else differs.

=head1 CONSTRUCTORS

Several functions take a string or numeric representation and generate
a corresponding date object.  The most general is C<new>, whose
argument list may be empty (returning the current date), a string in
format YYYY-MM-DD or YYYYMMDD, a list or arrayref of year, month, and
day number, or an existing date object.

=over 4

=item Date::Simple->new ([ARG, ...])

=item date ([ARG, ...])

    my $date = Date::Simple->new('1972-01-17');

The C<new> method will return a date object if the values passed in
specify a valid date.  (See above.)  If an invalid date is passed, the
method returns undef.  If the argument is invalid in form as opposed
to numeric range, C<new> dies.

The C<date> function provides the same functionality but must be
imported or qualified as C<Date::Simple::date>.  (To import all public
functions, do C<use Date::Simple (':all');>.)  This function returns
undef on all invalid input, rather than dying in some cases like
C<new>.

=item date_fmt (FMT,[ARG, ...])

Equivelent to C<date> but creates a Date::Simple::Fmt object instead. The
format is expected to be a valid POSIX::strftime format string.

=item date_iso ([ARG, ...])

Identical to C<date> but creates a Date::Simple::ISO object instead.

=item date_d8 ([ARG, ...])

Equivelent to C<date> but creates a Date::Simple::D8 object instead.

=item today()

Returns the current date according to C<localtime>.

B<Caution:> To get tomorrow's date (or any fixed offset from today),
do not use C<today + 1>.  Perl parses this as C<today(+1)>.  You need
to put empty parentheses after the function: C<today() + 1>.

=item ymd (YEAR, MONTH, DAY)

Returns a date object with the given year, month, and day numbers.  If
the arguments do not specify a valid date, undef is returned.

Example:

    use Date::Simple ('ymd');
    $pbd = ymd(1987, 12, 18);

=item d8 (STRING)

Parses STRING as "YYYYMMDD" and returns the corresponding date object,
or undef if STRING has the wrong format or specifies an invalid date.

Example:

    use Date::Simple ('d8');
    $doi = d8('17760704');

Mnemonic: The string matches C</\d{8}/>.  Also, "d8" spells "date", if
8 is expanded phonetically.

=back

=head1 INSTANCE METHODS

=over 4

=item DATE->next

    my $tomorrow = $today->next;

Returns an object representing tomorrow.

=item DATE->prev

   my $yesterday = $today->prev;

Returns an object representing yesterday.

=item DATE->year

    my $year  = $date->year;

Return the year of DATE as an integer.

=item DATE->month

    my $month = $date->month;

Return the month of DATE as an integer from 1 to 12.

=item DATE->day

    my $day   = $date->day;

Return the DATE's day of the month as an integer from 1 to 31.

=item DATE->day_of_week

Return a number representing DATE's day of the week from 0 to 6, where
0 means Sunday.

=item DATE->as_ymd

    my ($year, $month, $day) = $date->as_ymd;

Returns a list of three numbers: year, month, and day.

=item DATE->as_d8

Returns the "d8" representation (see C<d8>), like
C<$date-E<gt>format("%Y%m%d")>.

=item DATE->as_iso

Returns the ISO 8601 representation of the date (eg '2004-01-01'),
like C<$date-E<gt>format("%Y-%m-%d")>. This is in fact the default
overloaded stringification mechanism and is provided mostly so
other subclasses with different overloading can still do fast
ISO style date output.

=item DATE->as_str ([STRING])

=item DATE->format ([STRING])

=item DATE->strftime ([STRING])

These functions are equivalent.  Return a string representing the
date, in the format specified.  If you don't pass a parameter, the default
date format for the object is used if one has been specified, otherwise
uses the default date format for the class the object is a member of, or as
a last fallback uses the $Date::Simple::Standard_Format which is changeable,
but probably shouldn't be modified. See C<default_format> for details.

    my $change_date = $date->format("%d %b %y");
    my $iso_date1 = $date->format("%Y-%m-%d");
    my $iso_date2 = $date->format;

The formatting parameter is similar to one you would pass to
strftime(3).  This is because we actually do pass it to strftime to
format the date.  This may result in differing behavior across
platforms and locales and may not even work everywhere.

=item DATE->default_format ([FORMAT])

This method sets or gets the default_format for the DATE object or class
that it is called on.

=back

=head1 OPERATORS

Some operators can be used with Date::Simple instances.  If one side
of an expression is a date object, and the operator expects two date
objects, the other side is interpreted as C<date(ARG)>, so an array
reference or ISO 8601 string will work.

=over 4

=item DATE + NUMBER

=item DATE - NUMBER

You can construct a new date offset by a number of days using the C<+>
and C<-> operators.

=item DATE1 - DATE2

You can subtract two dates to find the number of days between them.

=item DATE1 == DATE2

=item DATE1 < DATE2

=item DATE1 <=> DATE2

=item DATE1 cmp DATE2

=item etc.

You can compare two dates using the arithmetic or string comparison
operators.  Equality tests (C<==> and C<eq>) return false when one of
the expressions can not be converted to a date.  Other comparison
tests die in such cases.  This is intentional, because in a sense, all
non-dates are not "equal" to all dates, but in no sense are they
"greater" or "less" than dates.

=item DATE += NUMBER

=item DATE -= NUMBER

You can increment or decrement a date by a number of days using the +=
and -= operators.  This actually generates a new date object and is
equivalent to C<$date = $date + $number>.

=item "$date"

You can interpolate a date instance directly into a string, in the
format specified by ISO 8601 (eg: 2000-01-17) for Date::Simple and
Date::Simple::ISO, for Date::Simple::D8 this is the same as calling
as_d8() on the object, and for Date::Simple::Fmt this is the same as
calling format() on the object.

=back

=head1 UTILITIES

=over 4

=item leap_year (YEAR)

Returns true if YEAR is a leap year.

=item days_in_month (YEAR, MONTH)

Returns the number of days in MONTH, YEAR.

=back

=over 4

=item leap_year (YEAR)

Returns true if YEAR is a leap year.

=item days_in_month (YEAR, MONTH)

Returns the number of days in MONTH, YEAR.

=back


=head1 AUTHOR

    Marty Pauley <marty@kasei.com>
    John Tobey <jtobey@john-edwin-tobey.org>
    Yves Orton <demerphq@hotmail.com>

=head1 COPYRIGHT

      Copyright (C) 2001  Kasei.
      Copyright (C) 2001,2002 John Tobey.
      Copyright (C) 2004 Yves Orton.

      This program is free software; you can redistribute it and/or
      modify it under the terms of either:

      a) the GNU General Public License;
         either version 2 of the License, or (at your option) any later
         version.  You should have received a copy of the GNU General
         Public License along with this program; see the file COPYING.
         If not, write to the Free Software Foundation, Inc., 59
         Temple Place, Suite 330, Boston, MA 02111-1307 USA

      b) the Perl Artistic License.

      This program is distributed in the hope that it will be useful,
      but WITHOUT ANY WARRANTY; without even the implied warranty of
      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=head1 SEE ALSO

L<Date::Simple::Fmt> L<Date::Simple::ISO> L<Date::Simple::D8>
and of course L<perl>

=cut
