#!/usr/bin/env perl
use strict;
use warnings;
( $#ARGV < 0) and die "Please specify a neutrino mass";
my $o_nu=shift; 
my $om=0.3;
my $hub=0.7;
my $ns=1.0;
if($#ARGV >= 2){
        $om=shift;
        $hub=shift;
        $ns=shift;
}
my $halo=1;

my $o_cdm=$om-0.05- $o_nu;
my $ol = 1 - $om;
my $m_nu=$o_nu*6/13.*100;
my $paramfile="default-params.ini";
my $newparams="params-nu$o_nu";
my $root = "out/nu$m_nu";
if(!$halo){ 
        $root=$root."-lin";
        $newparams = $newparams."-lin";
}
if($om != 0.3){ 
        $root = $root."om$om";
        $newparams = $newparams."om$om";
}
if($hub != 0.7){ 
        $root = $root."h$hub";
        $newparams = $newparams."h$hub";
}
if($ns != 1.0){ 
        $root = $root."ns$ns";
        $newparams = $newparams."ns$ns";
}
$newparams .=".ini";
my @red = (99,49,9, 4,3,2,1,0.5,0.3, 0.2,0.1,0);
$hub*=100;
#Read in template parameter file
open(my $INHAND, "<","$paramfile") or die "Could not open $paramfile for reading!";
open(my $OUTHAND, ">","$newparams") or die "Could not open $newparams for writing!";
while(<$INHAND>){
        s/^\s*output_root\s*=\s*[\w\/\.-]*/output_root=$root/i;
        #Set neutrino mass
        s/^\s*use_physical\s*=\s*[\w\/.-]*/use_physical = F/i;
        s/^\s*omega_cdm\s*=\s*[\w\/.-]*/omega_cdm=$o_cdm/i;
        s/^\s*omega_baryon\s*=\s*[\w\/.-]*/omega_baryon=0.05/i;
        s/^\s*omega_lambda\s*=\s*[\w\/.-]*/omega_lambda=$ol/i;
        s/^\s*omega_neutrino\s*=\s*[\w\/.-]*/omega_neutrino=$o_nu/i;
        s/^\s*hubble\s*=\s*[\w\/.-]*/hubble = $hub/i;
        #Set things we always need here
        s/^\s*massless_neutrinos\s*=\s*[\w\/.-]*/massless_neutrinos = 0/i;
        s/^\s*massive_neutrinos\s*=\s*[\w\/.-]*/massive_neutrinos = 3.04/i;
        s/^\s*get_transfer\s*=\s*[\w\/.-]*/get_transfer = T/i;
        s/^\s*do_nonlinear\s*=\s*[\w\/.-]*/do_nonlinear = $halo/i;
        #Set initial conditions; scalar_amp is to give sigma_8 = 0.878 at z=0 with nu=0.
        #Pivot irrelevant as n_s = 1
        s/^\s*initial_power_num\s*=\s*[\w\/.-]*/initial_power_num = 1/i;
        s/^\s*scalar_amp\(1\)\s*=\s*[\w\/.-]*/scalar_amp(1) = 2.27e-9/i;
        s/^\s*scalar_spectral_index\(1\)\s*=\s*[\w\/.-]*/scalar_spectral_index(1) = $ns/i;
        s/^\s*scalar_nrun\(1\)\s*=\s*[\w\/.-]*/scalar_nrun(1) = 0/i;
        #Set up output
        s/^\s*transfer_kmax\s*=\s*[\w\/.-]*/transfer_kmax = 30/i;
        my $nout=$#red+1;
        s/^\s*transfer_num_redshifts\s*=\s*[\w\/.-]*/transfer_num_redshifts = $nout/i;
        s/^\s*transfer_interp_matterpower\s*=\s*[\w\/.-]*/transfer_interp_matterpower = T/i;
        #Output files set later.
        for(my $i=0; $i<=$#red; $i++){
                s/^\s*transfer_redshift\($i\)\s*=\s*[\w\/.-]*//i;
                s/^\s*transfer_filename\($i\)\s*=\s*[\w\/.-]*//i;
                s/^\s*transfer_matterpower\($i\)\s*=\s*[\w\/.-]*//i;
        }
        #Write to new file.
        print $OUTHAND $_;
}
        #Set output files
        print $OUTHAND "\n#Transfer output files\n";
        for(my $i=0; $i<=$#red; $i++){
                print $OUTHAND "transfer_redshift(".($i+1).") = $red[$i]\n";
                print $OUTHAND "transfer_filename(".($i+1).") = transfer_$red[$i].dat\n";
                print $OUTHAND "transfer_matterpower(".($i+1).") = matterpow_$red[$i].dat\n";
        }
close($INHAND);
close($OUTHAND);

