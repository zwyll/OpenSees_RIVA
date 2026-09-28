# Run with OpenSees script.tcl /absolute/path/columnprobe.dylib /output/directory
load [lindex $argv 0]
set out [lindex $argv 1]
file mkdir $out
proc cube {latch} {
    wipe
    model BasicBuilder -ndm 3 -ndf 4
    set command {nDMaterial RIVASAND02 1 0.90 900 1.25 1.125 122.44207260468994 0.945 0.025 0.78 0.51 4 4 0.65 -rho 2.073519843851659 -nSub 20 -stage 1 -initialStress -100 -100 -100 0 0 0 -noBiasVolume}
    if {$latch} { lappend command -reversalLatch }
    {*}$command
    foreach {n x y z} {1 0 0 0 2 1 0 0 3 1 1 0 4 0 1 0 5 0 0 1 6 1 0 1 7 1 1 1 8 0 1 1} {
        node $n $x $y $z
    }
    element SSPbrickUP 1 1 2 3 4 5 6 7 8 1 550000 1 1e-5 1e-5 1e-5 0.537 2.381e-6 0 0 0
    columnProbe update
    columnProbe commit
}
proc shear {gamma} {
    foreach n {5 6 7 8} { setNodeDisp $n 1 $gamma }
    columnProbe update
}
set f [open "$out/trial_path.csv" w]
puts $f "latch,speculative_loading_trial,final_shear_strain,tau_xz_kPa,reversals"
foreach latch {0 1} {
    cube $latch
    for {set i 1} {$i<=100} {incr i} {
        shear [expr {$i*0.0001}]
        columnProbe commit
    }
    foreach speculative {0 1} {
        columnProbe revert
        if {$speculative} { shear 0.0101 }
        shear 0.0099
        puts $f "$latch,$speculative,[lindex [eleResponse 1 strain] 5],[lindex [eleResponse 1 stress] 5],[eleResponse 1 reversals]"
    }
}
close $f

cube 0
rayleigh 0 0 0 1
columnProbe commit
set initial [columnProbe matrix 1 tangent]
foreach n {5 6 7 8} { setNodeDisp $n 3 -0.001 }
columnProbe update
set trial [columnProbe matrix 1 tangent]
set damping [columnProbe matrix 1 damping]
rayleigh 0 0 0 0
set hydraulic [columnProbe matrix 1 damping]
set f [open "$out/damping_contract.csv" w]
puts $f "row,col,committed_stiffness,trial_stiffness,requested_committed_damping,no_rayleigh_damping"
set errorSquared 0.0
set referenceSquared 0.0
set committedErrorSquared 0.0
set committedReferenceSquared 0.0
set hydraulicError 0.0
for {set i 0} {$i<32} {incr i} {
    for {set j 0} {$j<32} {incr j} {
        set k [expr {32*$i+$j}]
        puts $f "$i,$j,[lindex $initial $k],[lindex $trial $k],[lindex $damping $k],[lindex $hydraulic $k]"
        if {$i%4!=3 && $j%4!=3} {
            set errorSquared [expr {$errorSquared+pow([lindex $damping $k]-[lindex $trial $k],2)}]
            set referenceSquared [expr {$referenceSquared+pow([lindex $trial $k],2)}]
            set committedErrorSquared [expr {$committedErrorSquared+pow([lindex $damping $k]-[lindex $initial $k],2)}]
            set committedReferenceSquared [expr {$committedReferenceSquared+pow([lindex $initial $k],2)}]
        } else {
            set hydraulicError [expr {max($hydraulicError,abs([lindex $damping $k]-[lindex $hydraulic $k]))}]
        }
    }
}
close $f
set dampingRelativeError [expr {sqrt($errorSquared/$referenceSquared)}]
set committedRelativeError [expr {sqrt($committedErrorSquared/$committedReferenceSquared)}]
cube 0
setNodeDisp 5 1 Inf
set rejected [catch {columnProbe update}]
set f [open "$out/checks.csv" w]
puts $f "trial_damping_relative_error,committed_damping_relative_error,hydraulic_damping_error,invalid_material_trial_rejected"
puts $f "$dampingRelativeError,$committedRelativeError,$hydraulicError,$rejected"
close $f
if {[lindex $argv 2] eq "strict"} {
    if {$dampingRelativeError>1.e-12} { error "Restored trial-stiffness damping contract failed" }
    if {$committedRelativeError<1.0} { error "Probe did not distinguish trial from committed stiffness" }
    if {$hydraulicError!=0.0} { error "Rayleigh term changed fluid/coupling blocks" }
    if {!$rejected} { error "Element silently accepted a rejected material trial" }
}
puts "CONTRACT_PROBES_COMPLETE $out"
wipe
