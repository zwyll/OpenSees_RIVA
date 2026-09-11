# Reuse the existing gravity/admission/stage tests with only TanType changed.
proc continuumNDMaterial {name args} {
    if {$name in {RIVASand RIVASAND02}} {lappend args -TanType 1}
    uplevel 1 [list nDMaterial $name {*}$args]
}
proc sourceContinuumTest {path} {
    set fp [open $path r]; set script [read $fp]; close $fp
    # BasicBuilder creates/recreates the native command at each model/wipe;
    # dispatch only these material calls through our opt-in test wrapper.
    uplevel #0 [string map {"nDMaterial " "continuumNDMaterial "} $script]
}
set here [file dirname [info script]]
sourceContinuumTest [file join $here RIVASand_stage_activation.tcl]
sourceContinuumTest [file join $here .. .. RIVASAND02 tests RIVASAND02_stage_activation.tcl]
if {[lindex [eleResponse 3 TanType] 0]!=1} {error "stage-2 material lost tangent selector"}
puts "PASS: continuum selector through elastic gravity, admission, and nonlinear/dynamic stages"
