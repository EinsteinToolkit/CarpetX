#!/usr/bin/perl

use strict;
use warnings;

my %methods = (
  "constant" => {"name" => "constant", 
                 "order" => 0},
  "euler" => {"name" => "Euler",
              "order" => 1},
  "rk2" => {"name" => "rk2",
            "order" => 2},
  "ssprk3" => {"name" => "ssprk3",
               "order" => 3},
  "rk4" => {"name" => "rk4",
            "order" => 4},
);

for my $method (keys %methods) {
  my $lines = <<EOF;
ActiveThorns = "
    CarpetX
    IOUtil
    ODESolvers
    TestODESolvers
"

Cactus::presync_mode = "mixed-error"

CarpetX::ncells_x = 1
CarpetX::ncells_y = 1
CarpetX::ncells_z = 1

CarpetX::blocking_factor_x = 1
CarpetX::blocking_factor_y = 1
CarpetX::blocking_factor_z = 1

CarpetX::ghost_size = 0

CarpetX::dtfac = 0.001
Cactus::cctk_itlast = 10

ODESolvers::method = "$methods{$method}->{name}"

# test something not exactly solvable by method
TestODESolvers::order = $methods{$method}->{order} + 1

IO::out_dir = \$parfile
IO::out_fileinfo = "axis labels"
IO::parfile_write = "no"

IO::out_every = 1
CarpetX::out_metadata = no
CarpetX::out_norm_omit_unstable = yes
CarpetX::out_norm_vars = "all"
EOF
  my $fn = "test-$method.par";
  open (my $fh, ">", $fn) or die "Could not open '$fn': $!";
  print $fh $lines;
  close $fh or die "Could not write to '$fn': $!";
}
