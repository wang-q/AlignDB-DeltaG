# NAME

AlignDB::DeltaG - Calculate deltaG of polymer DNA sequences

# SYNOPSIS

- Normal use

        use AlignDB::DeltaG
        my $obj = AlignDB::DeltaG->new(
            temperature => 37,
            salt_conc   => 1,
        );
        my $seq = "TAACAAGCAATGAGATAGAGAAAGAAATATATCCA";
        print "$seq deltaG: ", $obj->polymer_deltaG($seq), "\n";

- Reset conditionss

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

# DESCRIPTION

`AlignDB::DeltaG` is a simple class to calculate deltaG of polymer DNA sequences using the NN model.

In the near future, it may be extanded to calculate oligonucleotide thermodynamics.

## Reference

    1. SantaLucia J, Jr. 2004. Annu Rev Biophys Biomol Struct;
    2. SantaLucia J, Jr. 1998. Proc Natl Acad Sci U S A;

# ATTRIBUTES

## temperature

get/set temperature, Default: 37.0 degree centigrade

## salt\_conc

salt concentration, Default: 1 \[Na+\], in M
should be above 0.05 M and below 1.1 M

## deltaH

enthalpy, isa HashRef

## deltaS

entropy (cal/K.mol), isa HashRef

## deltaG

free energy, isa HashRef

# METHODS

## BUILD

rebuild the object by the new temperature and/or salt\_conc values

## polymer\_deltaG

my $dG = $obj->polymer\_deltaG($seq);
Calculate deltaG of a given sequence
This method is the main calculating sub.

## \_load\_thermodynamic\_data

my ($deltaH, $deltaS) = $self->\_load\_thermodynamic\_data;
load thermodynamic data come from references

## \_init\_deltaG

my $deltaG = $self->\_init\_deltaG;
recalculate deltaG by the new temperature and salt\_conc values

## \_rev\_com

my $rc\_polymer = $self->\_rev\_com($polymer);
Reverse and complementary a sequence

# AUTHOR

Qiang Wang &lt;wang-q@outlook.com>

# COPYRIGHT AND LICENSE

This software is copyright (c) 2008 by Qiang Wang.

This is free software; you can redistribute it and/or modify it under
the same terms as the Perl 5 programming language system itself.
