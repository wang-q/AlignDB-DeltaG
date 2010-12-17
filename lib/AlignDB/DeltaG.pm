package AlignDB::DeltaG;
# ABSTRACT: Calculate deltaG of polymer DNA sequences

use Moose;
use Carp;

=attr temperature

get/set temperature, Default: 37.0 degree centigrade

=cut

has 'temperature' => ( is => 'rw', isa => 'Num', default => 37.0, );

=attr salt_conc

salt concentration, Default: 1 [Na+], in M
should be above 0.05 M and below 1.1 M 
   
=cut
   
has 'salt_conc'   => ( is => 'rw', isa => 'Num', default => 1.0, );

=attr deltaH

enthalpy, isa HashRef

=cut

has 'deltaH' => ( is => 'ro', isa => 'HashRef', );

=attr deltaS

entropy (cal/K.mol), isa HashRef

=cut

has 'deltaS' => ( is => 'ro', isa => 'HashRef', );

=attr deltaG

free energy, isa HashRef

=cut

has 'deltaG' => ( is => 'ro', isa => 'HashRef', );

=method BUILD

      Usage : $obj->BUILD
    Purpose : rebuild the object by the new temperature and/or salt_conc values
    Returns : None
 Parameters : None
     Throws : no exceptions
   Comments : The BUILD method is called by Moose::Object::BUILDALL, which is
            : called by Moose::Object::new. So it is also the constructor
            : method. 
   See Also : n/a

=cut

sub BUILD {
    my $self = shift;

    # Load thermodynamic data
    my ( $deltaH, $deltaS ) = $self->_load_thermodynamic_data;
    $self->{deltaH} = $deltaH;
    $self->{deltaS} = $deltaS;

    # Recalculate the deltaG hash on current temperature and salt conditions
    my $deltaG = $self->_init_deltaG;
    $self->{deltaG} = $deltaG;

    return;
}

=method polymer_deltaG

      Usage : my $dG = $obj->polymer_deltaG($seq);
    Purpose : Calculate deltaG of a given sequence
    Returns : Num
 Parameters : Str
     Throws : no exceptions
   Comments : This method is the main calculating one.
   See Also : n/a

=cut

sub polymer_deltaG {
    my ( $self, $polymer ) = @_;

    $polymer = uc $polymer;
    return if $polymer =~ /[^AGCT]/;

    my $deltaG = $self->deltaG;

    my $polymer_len = length $polymer;
    my $dG          = 0;

    # calculate deltaG
    foreach ( 0 .. $polymer_len - 2 ) {
        my $nn = substr( $polymer, $_, 2 );
        $dG += $deltaG->{$nn};
    }

    # terminal correction
    my $init_terminal = "init" . substr( $polymer, 0, 1 );
    $dG += $deltaG->{$init_terminal};

    my $end_terminal = "init" . substr( $polymer, -1, 1 );
    $dG += $deltaG->{$end_terminal};

    # Symmetry correction
    my $rc_polymer = $self->_rev_com($polymer);
    if ( $polymer eq $rc_polymer ) {
        $dG += $deltaG->{sym};
    }

    return $dG;
}

=method _load_thermodynamic_data

      Usage : my ($deltaH, $deltaS) = $self->_load_thermodynamic_data;
    Purpose : load thermodynamic data come from references
    Returns : (HashRef, HashRef)
 Parameters : None
     Throws : no exceptions
   Comments : None
   See Also : n/a

=cut

sub _load_thermodynamic_data {
    my ($self) = @_;

    #-------------------#
    # deltaH (kcal/mol)
    #-------------------#
    my %deltaH = qw{
        AA -7.6 TT -7.6
        AT -7.2
        TA -7.2
        CA -8.5 TG -8.5
        GT -8.4 AC -8.4
        CT -7.8 AG -7.8
        GA -8.2 TC -8.2
        CG -10.6
        GC -9.8
        GG -8.0 CC -8.0
        initC 0.2 initG 0.2
        initA 2.2 initT 2.2
        sym 0.0
    };

    #--------------------#
    # deltaS (cal/K.mol)
    #--------------------#
    my %deltaS = qw{
        AA -21.3 TT -21.3
        AT -20.4
        TA -21.3
        CA -22.7 TG -22.7
        GT -22.4 AC -22.4
        CT -21.0 AG -21.0
        GA -22.2 TC -22.2
        CG -27.2
        GC -24.4
        GG -19.9 CC -19.9
        initC -5.7 initG -5.7
        initA 6.9 initT 6.9
        sym -1.4
    };

    return ( \%deltaH, \%deltaS );
}

=method _init_deltaG

      Usage : my $deltaG = $self->_init_deltaG;
    Purpose : recalculate deltaG by the new temperature and salt_conc values
    Returns : HashRef
 Parameters : None
     Throws : no exceptions
   Comments : None
   See Also : n/a

=cut

sub _init_deltaG {
    my ($self) = @_;

    # dG = dH - TdS, and dS is dependent on the salt concentration
    my $temperature = $self->temperature;
    my $salt_conc   = $self->salt_conc;
    my $deltaH      = $self->deltaH;
    my $deltaS      = $self->deltaS;

    my %deltaG = qw{
        initC 1.96
        initG 1.96
        initA 0.05
        initT 0.05
        sym 0.43
    };

    # the length of each NN dimer is 2, therefore the modifier is 1
    # total sodium concentration should be above 0.05 M and below 1.1 M
    my $entropy_adjust = ( 0.368 * log($salt_conc) );

    foreach my $key ( keys %{$deltaH} ) {

        # the length of each monomer is 1, thus the modifier of dS is 0
        # and the values are precalulated
        next if $key =~ /init|sym/;

        my $dS = $deltaS->{$key} + $entropy_adjust;
        my $dG = $deltaH->{$key}
            - ( ( 273.15 + $temperature ) * ( $dS / 1000 ) );
        $deltaG{$key} = $dG;
    }

    return \%deltaG;
}

=method _rev_com

      Usage : my $rc_polymer = $self->_rev_com($polymer);
    Purpose : Reverse and complementary a sequence
    Returns : Str
 Parameters : Str
     Throws : no exceptions
   Comments : None
   See Also : n/a

=cut

sub _rev_com {
    my ( $self, $sequence ) = @_;

    $sequence = reverse $sequence;    # reverse
    $sequence =~ tr/ACGTMRWSYKVHDBN/TGCAKYSWRMBDHVN/;    # complement

    return $sequence;
}

1;    # Magic true value required at end of module

__END__

=head1 SYNOPSIS

=over 2

=item Normal use

    use AlignDB::DeltaG
    my $obj = AlignDB::DeltaG->new(
        temperature => 37,
        salt_conc   => 1,
    );
    my $seq = "TAACAAGCAATGAGATAGAGAAAGAAATATATCCA";
    print "$seq deltaG: ", $obj->polymer_deltaG($seq), "\n";

=item Reset conditionss

    use AlignDB::DeltaG;
    # default value:
    #   temperature => 37,
    #   salt_conc   => 1,
    my $obj = AlignDB::DeltaG->new;
    $obj->temperature(30);
    $obj->salt_conc(0.1);
    $obj->BUILD;
    my $seq = "TAACAAGCAATGAGATAGAGAAAGAAATATATCCA";
    print "$seq deltaG: ", $obj->polymer_deltaG($seq), "\n";

=back

=head1 DESCRIPTION

C<AlignDB::DeltaG> is a simple class to calculate deltaG of polymer DNA
sequences using the NN model.

In the near future, it may be extanded to calculate oligonucleotide
thermodynamics.

=head2 Reference

 1. SantaLucia J, Jr. 2004. Annu Rev Biophys Biomol Struct; 
 2. SantaLucia J, Jr. 1998. Proc Natl Acad Sci U S A;

=cut
