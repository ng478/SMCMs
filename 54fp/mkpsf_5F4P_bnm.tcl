# VMD/psfgen script for generating psf/pdb pair for PDB 54FP
#  monomeric HIV-1 gp120 bound to bnm-III-170
#
# natasha gupta (c) 2018
# drexel university
# chemical and biological engineering
#
# check for base directory name variable;
# if not set, use default
if {![info exists PSFGEN_BASEDIR]} {
  # see if user set an environment variable
  if {[info exists env(PSFGEN_BASEDIR)]} {
      set PSFGEN_BASEDIR $env(PSFGEN_BASEDIR)
  } else {
      set PSFGEN_BASEDIR $env(HOME)/research/psfgen
  }
}
# load some custom TcL procedures to set coordinates correctly
source ${PSFGEN_BASEDIR}/src/loopmc.tcl
set LOCALFILES {}


mol new 5f4p.pdb
set molid [molinfo top get id]

# save heavy-atom coordinates for the four continuous potein fragments
set p1 [atomselect top "chain A and protein and resid 49 to 300"]
$p1 writepdb p1.pdb; lappend LOCALFILES p1.pdb
set p2 [atomselect top "chain A and protein and resid 325 to 457"]
$p2 writepdb p2.pdb; lappend LOCALFILES p2.pdb
set p3 [atomselect top "chain A and protein and resid 464 to 490"]
$p3 writepdb p3.pdb; lappend LOCALFILES p3.pdb

#
set n317pos [cacoIn_nOut 300 A 0]
set n458pos [cacoIn_nOut 457 A 0]
set n491pos [cacoIn_nOut 490 A 0]

# save the waters that are assigned to chain A that are close to the protein or ligand
set w [atomselect top "chain A and water and same residue as within 4.0 of (protein and chain A)"]
$w writepdb wA.pdb
lappend LOCALFILES wA.pdb

# save glycans
#set g [atomselect top "chain A and not water and not protein"]
#[atomselect top "chain A and resname NAG"] set resname BGNA
#[atomselect top "chain A and resname MAN"] set resname AMAN
#[atomselect top "chain A and resname BMA"] set resname BMAN
#[atomselect top "chain A and resname FUC"] set resname AFUC
#[atomselect top "chain A and resname GAL"] set resname BGAL
#$g writepdb gA.pdb
#lappend LOCALFILES gA.pdb

mol delete $molid


# Now we begin to build the composite PSF and PDB

package require psfgen

# Use the requisite charmm36 (july 2016) topology files
topology $env(HOME)/charmm/toppar/top_all36_prot.rtf
topology $env(HOME)/charmm/toppar/top_all36_carb_namd_cfa.rtf
topology $env(HOME)/charmm/toppar/stream/carb/toppar_all36_carb_glycopeptide.str
topology $env(HOME)/charmm/toppar/toppar_water_ions_namd.str

pdbalias atom ILE CD1 CD

# This alias assumes ALL histidines in the protein 
# are deprotonated at the delta position
pdbalias residue HIS HSD
pdbalias residue NAG BGLC
pdbalias atom BGLC C7 C
pdbalias atom BGLC O7 O
pdbalias atom BGLC C8 CT

# now we build the protein, including some
# missing residues
segment A {
  pdb p1.pdb
  residue 317 ASN
  residue 318 GLY
  residue 319 GLY
  residue 320 SER
  residue 321 GLY
  residue 322 SER
  residue 323 GLY
  residue 324 GLY
  pdb p2.pdb
  residue 458 GLY
  residue 459 GLY
  residue 460 ASN
  residue 461 ASP
  residue 462 ASP
  residue 463 ASN
  pdb p3.pdb
  residue 491 ILE
  
}

coordpdb p1.pdb A
coordpdb p2.pdb A
coordpdb p3.pdb A

# properly set the cartesian coordinates
# of the N's of the N-terminal
# missing residues of each missing segment
coord A 317 N $n317pos
coord A 458 N $n458pos
coord A 491 N $n491pos

# now we build the glycan segments
#segment GA {
#  pdb gA.pdb
#}
#coordpdb gA.pdb GA

# now we build the crystal water segments
pdbalias residue HOH TIP3
pdbalias atom TIP3 O OH2
segment WTXA {
  auto none
  pdb wA.pdb
}
coordpdb wA.pdb WTXA

# Now we create the disulfide bonds
patch DISU A:54 A:74
patch DISU A:119 A:205
patch DISU A:218 A:247
patch DISU A:228 A:239
patch DISU A:296 A:331
patch DISU A:378 A:445
patch DISU A:385 A:418

# Add glycan linkages
# patch NGLB A:241 GA:501
# patch NGLB A:262 GA:502
# patch NGLB A:276 GA:503
# patch NGLB A:289 GA:504
# patch NGLB A:295 GA:505
# patch NGLB A:356 GA:506
# patch NGLB A:386 GA:507
# patch NGLB A:394 GA:508
# patch NGLB A:448 GA:509

guesscoord

regenerate angles dihedrals

writepsf "my_5f4p.psf"
writepdb "my_5f4p_rawloops.pdb"
resetpsf

lappend LOCALFILES my_54fp_rawloops.pdb

mol new my_5f4p.psf
mol addfile my_5f4p_rawloops.pdb
set a [atomselect top all]
set molid [molinfo top get id]

set nc 1000
set rcut 3.0
set temperature 2.5
set k 10.0
set r0 1.5
set bg [atomselect ${molid} "noh"]                                                                                                                                   
set residueList [[atomselect ${molid} "chain A and resid 317 to 324 and name CA"] get residue]
do_loop_mc ${residueList} A ${molid} ${k} ${r0} ${bg} ${rcut} ${nc} ${temperature} [irand_dom 1000 9999]
set residueList [[atomselect ${molid} "chain A and resid 458 to 463 and name CA"] get residue]
do_loop_mc ${residueList} A ${molid} ${k} ${r0} ${bg} ${rcut} ${nc} ${temperature} [irand_dom 1000 9999]
set residueList [[atomselect ${molid} "chain A and resid 491 and name CA"] get residue]
do_loop_mc ${residueList} A ${molid} ${k} ${r0} ${bg} ${rcut} ${nc} ${temperature} [irand_dom 1000 9999]

$a writepdb "my_54fp_mcOut.pdb"

set fix [atomselect top "protein and noh and not (resid 300 to 325 457 to 464)"]
$a set beta 0
$fix set beta 1
set wat [atomselect top "name OH2"]
$wat set beta 1
set glyhv [atomselect top "not water and not protein and noh"]
$glyhv set beta 1
$a writepdb "my_54fp_fix.pdb"

# clean up
foreach f $LOCALFILES {
  exec rm $f
}

puts "Generated my_54fp.psf, my_3tgq_mc.pdb, and my_3tgq_fix.pdb."

exit
