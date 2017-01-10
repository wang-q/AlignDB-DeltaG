requires 'Moose';
requires 'YAML::Syck', '1.29';
requires 'perl', '5.010001';

on test => sub {
    requires 'Test::More', 0.88;
};
