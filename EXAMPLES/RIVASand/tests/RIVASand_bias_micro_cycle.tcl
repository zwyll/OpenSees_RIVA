# Regression test for the bias_reversible_volume amplitude/pressure gates.
#
# Two configurations that were pathological before the gates (see
# EXAMPLES/RIVASand/BIAS_VOLUME_LOW_CONFINEMENT.md):
#   case A: sigma_v0 = 25 kPa, static bias alpha = 0.09, tau_cyc = 0.005 kPa
#           (micro-amplitude cycling on a biased point)
#   case B: sigma_v0 = 2.5 kPa, alpha = 0.09, CSR = 0.05
#           (low-confinement shallow-slope state)
# Pre-gate behavior: sigma'_v oscillates with a GROWING envelope driven by the
# phase-keyed bias_reversible_volume target (case A: 25 -> 13.5/47.8 kPa in
# five cycles with shear strain frozen at ~3.5e-5). Post-gate: bounded.
#
# PASS criterion per case: max |sigma'_v - sigma'_v0| / sigma'_v0 < 0.15
# over 5 cycles of constant-volume stress-controlled DSS.
wipe

proc runCase {label sigv csr} {
    wipe
    model BasicBuilder -ndm 3 -ndf 3
    set Dr 0.641
    set alpha 0.09
    nDMaterial RIVASand 1 $Dr 1.25 1.125 122.44207260468994 0.945 0.025 \
        0.78 0.51 10.0 1.5 0.65 -nSub 8 -stressScale 1.0 -tangentPMin 2.0 -stage 0
    node 1 0.0 0.0 0.0
    node 2 1.0 0.0 0.0
    node 3 1.0 1.0 0.0
    node 4 0.0 1.0 0.0
    node 5 0.0 0.0 1.0
    node 6 1.0 0.0 1.0
    node 7 1.0 1.0 1.0
    node 8 0.0 1.0 1.0
    element SSPbrick 1 1 2 3 4 5 6 7 8 1 0.0 0.0 0.0
    fix 1 1 1 1
    fix 2 1 1 1
    fix 3 1 1 1
    fix 4 1 1 1
    foreach n {5 6 7 8} { fix $n 0 1 0 }
    equalDOF 5 6 1
    equalDOF 5 7 1
    equalDOF 5 8 1
    pattern Plain 1 Linear {
        foreach n {5 6 7 8} {
            eval "load $n [expr {$alpha*$sigv/4.0}] 0.0 [expr {-$sigv/4.0}]"
        }
    }
    constraints Transformation
    numberer Plain
    system ProfileSPD
    test NormDispIncr 1.0e-8 100 0
    algorithm KrylovNewton
    integrator LoadControl 0.02
    analysis Static
    if {[analyze 50] != 0} { puts "FAIL $label consolidation"; return 1 }
    loadConst -time 0.0
    updateMaterialStage -material 1 -stage 1
    integrator LoadControl 0.0
    test NormDispIncr 1.0e-7 100 0
    for {set i 1} {$i <= 10} {incr i} {
        if {[analyze 1] != 0} { puts "FAIL $label hold"; return 1 }
    }
    foreach n {5 6 7 8} { sp $n 3 [nodeDisp $n 3] }
    set tauAmp [expr {$csr*$sigv}]
    timeSeries Trig 2 0.0 5.0 1.0 -factor 1.0
    pattern Plain 2 2 {
        foreach n {5 6 7 8} { eval "load $n [expr {$tauAmp/4.0}] 0.0 0.0" }
    }
    integrator LoadControl 0.005
    set drift 0.0
    for {set i 1} {$i <= 1000} {incr i} {
        if {[analyze 1] != 0} { puts "FAIL $label cyclic step $i"; return 1 }
        set s [eleResponse 1 stress]
        set sv [expr {-[lindex $s 2]}]
        set d [expr {abs($sv - $sigv)/$sigv}]
        if {$d > $drift} { set drift $d }
    }
    if {$drift < 0.15} {
        puts "PASS $label: max |sigma_v - sigma_v0|/sigma_v0 = [format %.4f $drift]"
        return 0
    }
    puts "FAIL $label: sigma_v drift [format %.4f $drift] exceeds 0.15 (growing bias-volume oscillation)"
    return 1
}

set bad 0
incr bad [runCase "microAmp_sv25" 25.0 0.0002]
incr bad [runCase "lowConf_sv2.5" 2.5 0.05]
if {$bad == 0} { puts "PASS: RIVASand bias micro-cycle regression (2 cases)"; exit 0 }
puts "FAIL: $bad case(s)"
exit 1
