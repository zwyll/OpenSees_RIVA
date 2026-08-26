# Verify Tcl parsing and element-copy state for compressive, stress-preserving
# over-bound geostatic admission. Stress unit: kPa.

wipe
model BasicBuilder -ndm 3 -ndf 3

set matTag 8010
set Dr 0.641
set initialStress {-2.0 -2.0 -56.0 0.0 0.0 0.0}
nDMaterial RIVASand $matTag \
    $Dr 1.25 1.125 122.44207260468994 0.945 0.025 \
    0.78 0.51 10.0 1.5 0.65 \
    -nSub 1 -stressScale 1.0 -geostaticAdmission \
    -initialStress {*}$initialStress

node 1 0.0 0.0 0.0
node 2 1.0 0.0 0.0
node 3 1.0 1.0 0.0
node 4 0.0 1.0 0.0
node 5 0.0 0.0 1.0
node 6 1.0 0.0 1.0
node 7 1.0 1.0 1.0
node 8 0.0 1.0 1.0
element bbarBrick 1 1 2 3 4 5 6 7 8 $matTag

set stress [eleResponse 1 material 1 stress]
if {[llength $stress] != 6} {
    error "RIVASand admitted stress response is unavailable"
}
for {set i 0} {$i < 6} {incr i} {
    if {abs([lindex $stress $i]-[lindex $initialStress $i]) > 1.0e-12} {
        error "geostatic admission changed stress component $i: $stress"
    }
}

set state [eleResponse 1 material 1 state]
if {[llength $state] != 93} {
    error "RIVASand admitted state has [llength $state] values, expected 93"
}
set admissionRadius [lindex [eleResponse 1 material 1 geostaticAdmissionRadius] 0]
set admitted [lindex [eleResponse 1 material 1 geostaticAdmitted] 0]
if {$admitted != 1.0 || $admissionRadius <= 0.0} {
    error "over-bound state was not admitted: radius=$admissionRadius flag=$admitted"
}

set dbPath [file join [pwd] "RIVASand_geostatic_admission_[pid]"]
file delete -force $dbPath
database File $dbPath
save 1
restore 1
set restartedRadius [lindex \
    [eleResponse 1 material 1 geostaticAdmissionRadius] 0]
set restartedFlag [lindex [eleResponse 1 material 1 geostaticAdmitted] 0]
file delete -force $dbPath
if {$restartedFlag != 1.0 ||
    abs($restartedRadius-$admissionRadius) > 1.0e-12} {
    error "restart lost admission state: radius=$restartedRadius flag=$restartedFlag"
}

puts "PASS: Tcl geostatic admission preserves stress and restart state; radius=$admissionRadius"
