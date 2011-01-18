#!/usr/bin/env perl
use strict;
use warnings;
( $#ARGV < 1) and die "Please specify a neutrino mass";
my $o_nu=shift; 
my $halo=shift;
my $o_cdm=0.25- $o_nu;
my $m_nu=$o_nu*6/13.*100;
my $paramfile="default-params.ini";
my $newparams="params-nu$o_nu.ini";
if(!$halo){
$newparams="params-nu$o_nu-lin.ini";
}
my @red = (9, 4,3,2,1,0.5, 0.2,0);
#Read in template parameter file
open(my $INHAND, "<","$paramfile") or die "Could not open $paramfile for reading!";
open(my $OUTHAND, ">","$newparams") or die "Could not open $newparams for writing!";
while(<$INHAND>){
        if(!$halo){
        s/^\s*output_root\s*=\s*[\w\/\.-]*/output_root=out\/nu$m_nu-lin/i;
        }
        else{
        s/^\s*output_root\s*=\s*[\w\/\.-]*/output_root=out\/nu$m_nu/i;
        }
        #Set neutrino mass
        s/^\s*use_physical\s*=\s*[\w\/.-]*/use_physical = F/i;
        s/^\s*omega_cdm\s*=\s*[\w\/.-]*/omega_cdm=$o_cdm/i;
        s/^\s*omega_baryon\s*=\s*[\w\/.-]*/omega_baryon=0.05/i;
        s/^\s*omega_neutrino\s*=\s*[\w\/.-]*/omega_neutrino=$o_nu/i;
        s/^\s*hubble\s*=\s*[\w\/.-]*/hubble = 70/i;
        #Set things we always need here
        s/^\s*massless_neutrinos\s*=\s*[\w\/.-]*/massless_neutrinos = 0/i;
        s/^\s*massive_neutrinos\s*=\s*[\w\/.-]*/massive_neutrinos = 3.04/i;
        s/^\s*get_transfer\s*=\s*[\w\/.-]*/get_transfer = T/i;
        s/^\s*do_nonlinear\s*=\s*[\w\/.-]*/do_nonlinear = $halo/i;
        #Set initial conditions; scalar_amp is to give sigma_8 = 0.878 at z=0 with nu=0.
        #Pivot irrelevant as n_s = 1
        s/^\s*initial_power_num\s*=\s*[\w\/.-]*/initial_power_num = 1/i;
        s/^\s*scalar_amp\(1\)\s*=\s*[\w\/.-]*/scalar_amp(1) = 2.27e-9/i;
        s/^\s*scalar_spectral_index\(1\)\s*=\s*[\w\/.-]*/scalar_spectral_index(1) = 1.0/i;
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

