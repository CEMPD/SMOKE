package Date::Simple::ISO;
use Date::Simple 3;
use base qw/Date::Simple/;
use overload '""' => 'as_iso';    # sub { $_[0]->as_iso };

*EXPORT      = *Date::Simple::EXPORT;
*EXPORT_OK   = *Date::Simple::EXPORT_OK;
*EXPORT_TAGS = *Date::Simple::EXPORT_TAGS;

sub d8    { shift->_d8(@_); }
sub today { shift->_today(@_); }
sub ymd   { shift->_ymd(@_); }

1;

=head1 NAME

Date::Simple::ISO - Sub class of Date::Simple

=head1 SYNOPSIS

    use Date::Simple::ISO;

=head1 DESCRIPTION

This module is entirely identical to Date::Simple. It is included for completness
and self documenting sake.  IMO it is preferable to say

  my $obj=Date::Simple::ISO->new(...);

As this makes the implicit formatting of the object clear.

=head1 SEE ALSO

L<Date::Simple> for full documentation.

=cut
