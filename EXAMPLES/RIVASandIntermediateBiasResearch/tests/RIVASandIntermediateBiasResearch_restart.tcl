# Verify OpenSees sendSelf/recvSelf for the revision-3 research state.

proc configureAnalysis {} {
    wipeAnalysis
    constraints Transformation
    numberer RCM
    system BandGeneral
    test NormDispIncr 1.0e-12 12 0
    algorithm Newton
    integrator DisplacementControl 5 1 1.0e-12
    analysis Static
}

proc advanceCycle {firstStep lastStep pointsPerCycle amplitude} {
    set pi [expr {acos(-1.0)}]
    set previousGamma [expr {$amplitude*sin(2.0*$pi*double($firstStep-1)/$pointsPerCycle)}]
    for {set step $firstStep} {$step <= $lastStep} {incr step} {
        set gamma [expr {$amplitude*sin(2.0*$pi*double($step)/$pointsPerCycle)}]
        integrator DisplacementControl 5 1 [expr {$gamma-$previousGamma}]
        if {[analyze 1] != 0} {
            error "restart test failed at cyclic step $step"
        }
        set previousGamma $gamma
    }
}

proc assertVectorsClose {actual expected tolerance label} {
    if {[llength $actual] != [llength $expected]} {
        error "$label length [llength $actual] != [llength $expected]"
    }
    for {set i 0} {$i < [llength $actual]} {incr i} {
        set a [lindex $actual $i]
        set e [lindex $expected $i]
        if {abs($a-$e) > $tolerance} {
            error "$label differs at $i: $a versus $e"
        }
    }
}

wipe
model BasicBuilder -ndm 3 -ndf 3

set matTag 8010
set Dr [expr {(0.78-0.601)/(0.78-0.51)}]
nDMaterial RIVASandIntermediateBiasResearch $matTag \
    $Dr 1.25 1.125 122.44207260468994 0.945 0.025 \
    0.78 0.51 10.0 1.5 0.65 -nSub 4 -stressScale 1.0 \
    -reversalLatch \
    -stage 1 -initialStress -19.4 -19.4 -40.0 0.0 0.0 10.0

foreach {tag x y z} {
    1 0 0 0  2 1 0 0  3 1 1 0  4 0 1 0
    5 0 0 1  6 1 0 1  7 1 1 1  8 0 1 1
} {
    node $tag $x $y $z
}
element bbarBrick 1 1 2 3 4 5 6 7 8 $matTag
foreach tag {1 2 3 4} { fix $tag 1 1 1 }
foreach tag {5 6 7 8} { fix $tag 0 1 1 }
equalDOF 5 6 1
equalDOF 5 7 1
equalDOF 5 8 1
pattern Plain 1 Linear { load 5 1.0 0.0 0.0 }
configureAnalysis

set copiedLatch [lindex [eleResponse 1 material 1 reversalLatch] 0]
if {$copiedLatch != 1.0} {
    error "element material copy lost -reversalLatch"
}

set pointsPerCycle 32
set amplitude 0.005
advanceCycle 1 32 $pointsPerCycle $amplitude
set checkpointState [eleResponse 1 material 1 state]
if {[llength $checkpointState] != 138 || [lindex $checkpointState 137] <= 0.0} {
    error "checkpoint did not contain the revision-3 plastic-activity state"
}

set databasePrefix [file join [pwd] RIVASandIntermediateBiasResearch_restart_db]
foreach path [glob -nocomplain ${databasePrefix}*] { file delete -force $path }
database File $databasePrefix
save 101

advanceCycle 33 64 $pointsPerCycle $amplitude
set baselineState [eleResponse 1 material 1 state]
set baselineStress [eleResponse 1 material 1 stress]

restore 101
configureAnalysis
set restoredLatch [lindex [eleResponse 1 material 1 reversalLatch] 0]
if {$restoredLatch != 1.0} {
    error "sendSelf/recvSelf lost -reversalLatch"
}
advanceCycle 33 64 $pointsPerCycle $amplitude
set restartedState [eleResponse 1 material 1 state]
set restartedStress [eleResponse 1 material 1 stress]

assertVectorsClose $restartedState $baselineState 1.0e-13 "restart state"
assertVectorsClose $restartedStress $baselineStress 1.0e-13 "restart stress"

wipe
foreach path [glob -nocomplain ${databasePrefix}*] { file delete -force $path }
puts "PASS: revision-3 RIVASandIntermediateBiasResearch save/restore and serialized reversal-latch continuity"
