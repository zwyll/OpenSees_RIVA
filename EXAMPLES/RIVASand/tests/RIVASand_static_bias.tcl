# Static-bias activation test: reproduce the SSI stage-transition failure at a
# single element. Elastic (stage-0) K0-type compression to SIGV with a static
# shear bias TAU = BIAS*SIGV applied on top, commit, activate stage 1 (the
# kernel initializes from the committed biased state), then HOLD the loads for
# a few equilibrium steps. A sustainable state converges with small strains; an
# unsustainable one flows / fails to converge — the single-element analogue of
# the foundation-edge Gauss point in an SSI gravity field.
#
# env: RIVA_DR (0.50), RIVA_SIGV (kPa, 25), RIVA_BIAS (0.30),
#      RIVA_APEX (kPa, 0 = frozen kernel), RIVA_ZETA/H/KD/M (row knobs)
wipe
model BasicBuilder -ndm 3 -ndf 3

proc envv {name def} { return [expr {[info exists ::env($name)] ? $::env($name) : $def}] }
set Dr    [envv RIVA_DR 0.50]
set SIGV  [envv RIVA_SIGV 25.0]
set BIAS  [envv RIVA_BIAS 0.30]
set APEX  [envv RIVA_APEX 0.0]
set LATCH [envv RIVA_LATCH 0]
set zeta  [envv RIVA_ZETA 0.05]
set hh    [envv RIVA_H 380.0]
set kd    [envv RIVA_KD 1.125]
set MM    [envv RIVA_M 1.25]

set cmd [list nDMaterial RIVASand 1 $Dr $MM $kd $hh 0.945 $zeta \
    0.78 0.51 10.0 1.5 0.65 -rho 1.965 -nSub 20 -stressScale 1.0 \
    -tangentPMin 2.0 -stage 0]
if {$APEX > 0.0} { lappend cmd -pMin $APEX }
if {$LATCH != 0} { lappend cmd -reversalLatch }
if {[envv RIVA_ADMIT 0] != 0} { lappend cmd -admitOverbound }
if {[envv RIVA_NUMTAN 0] != 0} { lappend cmd -numericalTangent }
if {[envv RIVA_BETAF 0.0] > 0.0} { lappend cmd -betaFloor [envv RIVA_BETAF 0.0] }
if {[envv RIVA_RECENTER 0] != 0} { lappend cmd -recenterActivation }
if {[envv RIVA_BETARES 0.0] > 0.0} { lappend cmd -betaReserve [envv RIVA_BETARES 0.0] }
eval $cmd

node 1 0.0 0.0 0.0
node 2 1.0 0.0 0.0
node 3 1.0 1.0 0.0
node 4 0.0 1.0 0.0
node 5 0.0 0.0 1.0
node 6 1.0 0.0 1.0
node 7 1.0 1.0 1.0
node 8 0.0 1.0 1.0
element SSPbrick 1 1 2 3 4 5 6 7 8 1 0.0 0.0 0.0

# base fixed; top: x free (tied), y fixed, z free (tied) — DSS kinematics
fix 1 1 1 1
fix 2 1 1 1
fix 3 1 1 1
fix 4 1 1 1
foreach n {5 6 7 8} { fix $n 0 1 0 }
equalDOF 5 6 1 3
equalDOF 5 7 1 3
equalDOF 5 8 1 3

# stage 0 (elastic): vertical SIGV + shear BIAS*SIGV together, 10 steps
set FV [expr {-$SIGV/4.0}]
set FH [expr {$BIAS*$SIGV/4.0}]
pattern Plain 1 Linear {
    foreach n {5 6 7 8} {
        eval "load $n $FH 0.0 $FV"
    }
}
constraints Transformation
numberer RCM
system BandGeneral
test NormDispIncr 1.0e-10 50 0
algorithm Newton
integrator LoadControl 0.1
analysis Static
if {[analyze 10] != 0} { puts "RESULT elastic-load FAILED"; exit }
loadConst -time 0.0

set s0 [eleResponse 1 stress]
puts "committed stage-0 stress: $s0"

updateMaterialStage -material 1 -stage 1

# hold the same loads: 20 zero-increment equilibrium steps
integrator LoadControl 0.0
test NormDispIncr 1.0e-8 100 0
set ok 0
set gmax 0.0
for {set i 1} {$i <= 20} {incr i} {
    if {[analyze 1] != 0} { set ok -1; break }
    set eps [eleResponse 1 strain]
    set g [expr {abs([lindex $eps 5])}]
    if {$g > $gmax} { set gmax $g }
    if {$g > 0.10} { set ok -2; break }
}
if {$ok != 0} {
    puts "RESULT DIVERGE-HOLD Dr=$Dr sigv=$SIGV bias=$BIAS apex=$APEX gmax=[format %.4g $gmax]"
    wipe
    exit
}

# phase 3: perturb — small strain-controlled shear push from the biased state
# (the single-element analogue of the first post-activation solver increment)
# RIVA_TRACE=1: use algorithm Linear (no iteration) with tiny steps and print
# the raw tau(gamma)/p(gamma) path — diagnoses the material response itself.
set pushTot [expr {1.0e-3*[envv RIVA_PUSHDIR 1.0]}]
set nPush 20
set TRACE [envv RIVA_TRACE 0]
pattern Plain 2 Linear { load 5 1.0 0.0 0.0 }
integrator DisplacementControl 5 1 [expr {$pushTot/$nPush}]
test NormDispIncr 1.0e-8 100 0
if {$TRACE != 0} {
    algorithm Newton
    test FixedNumIter 1
    set nPush 40
    integrator DisplacementControl 5 1 [expr {$pushTot/$nPush}]
}
set ok 0
for {set i 1} {$i <= $nPush} {incr i} {
    if {[analyze 1] != 0} { set ok -1; break }
    if {$TRACE != 0} {
        set s [eleResponse 1 stress]
        set e [eleResponse 1 strain]
        set pp [expr {-([lindex $s 0]+[lindex $s 1]+[lindex $s 2])/3.0}]
        puts "TRACE g=[format %.3e [lindex $e 5]] tau=[format %.4g [lindex $s 5]] p=[format %.4g $pp]"
        flush stdout
    }
}
set s1 [eleResponse 1 stress]
set pE [expr {-([lindex $s1 0]+[lindex $s1 1]+[lindex $s1 2])/3.0}]
set tE [lindex $s1 5]
set eps [eleResponse 1 strain]
set gE [lindex $eps 5]
set tag STABLE
if {$ok != 0} { set tag DIVERGE-PUSH }
if {$pE < 0.05} { set tag COLLAPSE }
puts "RESULT $tag Dr=$Dr sigv=$SIGV bias=$BIAS apex=$APEX pEnd=[format %.3g $pE] tauEnd=[format %.3g $tE] gEnd=[format %.3g $gE]"
wipe
