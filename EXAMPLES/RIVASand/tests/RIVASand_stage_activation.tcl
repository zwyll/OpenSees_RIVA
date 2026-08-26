# Verify that updateMaterialStage reaches every RIVASand copy held by
# multiple OpenSees elements and initializes each from its own committed
# effective stress.

wipe
model BasicBuilder -ndm 3 -ndf 4

set matTag 8002
set Dr [expr {(0.78-0.601)/(0.78-0.51)}]
nDMaterial RIVASand $matTag \
    $Dr 1.25 1.125 122.44207260468994 0.945 0.025 \
    0.78 0.51 10.0 1.5 0.65 -nSub 1 -stressScale 1.0 \
    -tangentPMin 2.0

node 1 0.0 0.0 0.0
node 2 1.0 0.0 0.0
node 3 1.0 1.0 0.0
node 4 0.0 1.0 0.0
node 5 0.0 0.0 1.0
node 6 1.0 0.0 1.0
node 7 1.0 1.0 1.0
node 8 0.0 1.0 1.0
node 9  0.0 0.0 2.0
node 10 1.0 0.0 2.0
node 11 1.0 1.0 2.0
node 12 0.0 1.0 2.0
element SSPbrickUP 1 1 2 3 4 5 6 7 8 $matTag \
    2.2e6 1.0 1.0e-8 1.0e-8 1.0e-8 0.601 1.0e-5
element SSPbrickUP 2 5 6 7 8 9 10 11 12 $matTag \
    2.2e6 1.0 1.0e-8 1.0e-8 1.0e-8 0.601 1.0e-5

fix 1 1 1 1 1
fix 2 1 1 1 1
fix 3 1 1 1 1
fix 4 1 1 1 1
fix 5 0 1 0 1
fix 6 0 1 0 1
fix 7 0 1 0 1
fix 8 0 1 0 1
fix 9  0 1 0 1
fix 10 0 1 0 1
fix 11 0 1 0 1
fix 12 0 1 0 1
equalDOF 5 6 1 3
equalDOF 5 7 1 3
equalDOF 5 8 1 3
equalDOF 9 10 1 3
equalDOF 9 11 1 3
equalDOF 9 12 1 3

pattern Plain 1 Linear {
    load 9 0.0 0.0 -1.0 0.0
}

constraints Transformation
numberer RCM
system BandGeneral
test NormDispIncr 1.0e-12 12 0
algorithm Newton
integrator DisplacementControl 9 3 -0.002
analysis Static
if {[analyze 1] != 0} {
    error "RIVASand stage-0 compression failed"
}

loadConst -time 0.0
updateMaterialStage -material $matTag -stage 1

set pressureAnchors {}
foreach eleTag {1 2} {
    set initializedState [eleResponse $eleTag state]
    if {[llength $initializedState] != 94} {
        error "RIVASand element $eleTag state response has [llength $initializedState] values, expected 94"
    }
    set pressureAnchor [lindex $initializedState 49]
    if {$pressureAnchor <= 0.0} {
        error "RIVASand element $eleTag copy was not initialized at stage activation"
    }
    set materialStage [lindex [eleResponse $eleTag stage] 0]
    if {$materialStage != 1.0} {
        error "RIVASand element $eleTag reports stage=$materialStage"
    }
    set pressureFloor [lindex [eleResponse $eleTag pressureFloor] 0]
    set tangentFloor [lindex [eleResponse $eleTag tangentPressureFloor] 0]
    if {abs($pressureFloor-0.001) > 1.0e-12 ||
        abs($tangentFloor-2.0) > 1.0e-12} {
        error "RIVASand element $eleTag floors are pMin=$pressureFloor tangentPMin=$tangentFloor"
    }
    lappend pressureAnchors $pressureAnchor
}

pattern Plain 2 Linear {
    load 9 1.0 0.0 0.0 0.0
}
integrator DisplacementControl 9 1 0.0005
if {[analyze 1] != 0} {
    error "RIVASand stage-1 shear step failed"
}

foreach eleTag {1 2} {
    set stress [eleResponse $eleTag stress]
    if {[llength $stress] != 6} {
        error "RIVASand element $eleTag stress response is unavailable after activation"
    }
}
puts "PASS: every RIVASand copy activated; pAnchors=$pressureAnchors kPa"
