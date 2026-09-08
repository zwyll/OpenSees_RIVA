# Exercise the public selector through cyclic element histories and restart.

proc materialCommand {options {initialStress {}}} {
    set Dr [expr {(0.78-0.601)/(0.78-0.51)}]
    set command [list nDMaterial RIVASAND02 8030 \
        $Dr 483.48301127222084 1.25 1.125 122.44207260468994 0.945 0.025 \
        0.78 0.51 10.0 1.5 0.65 -nSub 4 -stressScale 1.0]
    if {[llength $initialStress] != 0} {
        lappend command -stage 1 -initialStress {*}$initialStress
    }
    return [concat $command $options]
}

proc assertClose {actual expected tolerance label} {
    if {[llength $actual] != [llength $expected]} {
        error "$label: response lengths differ"
    }
    foreach a $actual e $expected {
        # The negated comparison also rejects non-finite differences.
        if {![expr {abs($a-$e) <= $tolerance}]} {
            error "$label: $a differs from $e"
        }
    }
}

proc assertMode {mode} {
    foreach key {BiasVolume biasVolume} {
        assertClose [eleResponse 1 material 1 $key] [list $mode] 0.0 $key
    }
    assertClose [eleResponse 1 material 1 noBiasVolume] \
        [list [expr {$mode == 1}]] 0.0 "legacy noBiasVolume response"
    assertClose [eleResponse 1 material 1 fieldBiasVolume] \
        [list [expr {$mode == 2}]] 0.0 "legacy fieldBiasVolume response"
    assertClose [eleResponse 1 material 1 reversalLatch] {1} 0.0 "reversalLatch"
}

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

proc advanceHistory {first last} {
    set history {}
    set pi [expr {acos(-1.0)}]
    for {set step $first} {$step <= $last} {incr step} {
        set gamma [expr {0.005*sin(2.0*$pi*$step/32.0)}]
        set previous [expr {0.005*sin(2.0*$pi*($step-1)/32.0)}]
        integrator DisplacementControl 5 1 [expr {$gamma-$previous}]
        if {[analyze 1] != 0} { error "cyclic step $step failed" }
        set state [eleResponse 1 material 1 state]
        if {[llength $state] != 139} { error "incomplete material state" }
        lappend history {*}[eleResponse 1 material 1 stress] {*}$state
    }
    return $history
}

proc runHistory {options mode initialStress databasePrefix} {
    wipe
    model BasicBuilder -ndm 3 -ndf 3
    {*}[materialCommand [concat {-reversalLatch} $options] $initialStress]
    foreach {tag x y z} {
        1 0 0 0  2 1 0 0  3 1 1 0  4 0 1 0
        5 0 0 1  6 1 0 1  7 1 1 1  8 0 1 1
    } { node $tag $x $y $z }
    element bbarBrick 1 1 2 3 4 5 6 7 8 8030
    foreach tag {1 2 3 4} { fix $tag 1 1 1 }
    foreach tag {5 6 7 8} { fix $tag 0 1 1 }
    foreach tag {6 7 8} { equalDOF 5 $tag 1 }
    pattern Plain 1 Linear { load 5 1.0 0.0 0.0 }
    configureAnalysis
    assertMode $mode
    set history [advanceHistory 1 32]
    if {$databasePrefix ne ""} {
        database File $databasePrefix
        save 301
    }
    set tail [advanceHistory 33 64]
    if {$databasePrefix ne ""} {
        restore 301
        configureAnalysis
        assertMode $mode
        assertClose [advanceHistory 33 64] $tail 1.0e-13 "mode $mode restart"
    }
    wipe
    return [concat $history $tail]
}

# Malformed inputs and contradictory aliases must not silently select a law.
set badInputs [list {-BiasVolume} {-BiasVolume -1} {-BiasVolume 3} \
    {-BiasVolume 1.5} {-BiasVolume word} {-BiasVolume -nSub 4} \
    {-BiasVolume 0 -BiasVolume 1} {-BiasVolume 2 -BiasVolume 1} \
    {-BiasVolume 0 -noBiasVolume} {-noBiasVolume -BiasVolume 0} \
    {-BiasVolume 0 -fieldBiasVolume} {-fieldBiasVolume -BiasVolume 0} \
    {-BiasVolume 1 -fieldBiasVolume} {-fieldBiasVolume -BiasVolume 1} \
    {-BiasVolume 2 -noBiasVolume} {-noBiasVolume -BiasVolume 2} \
    {-noBiasVolume -fieldBiasVolume} {-fieldBiasVolume -noBiasVolume}]
foreach options $badInputs {
    wipe
    model BasicBuilder -ndm 3 -ndf 3
    if {![catch {{*}[materialCommand $options]}]} {
        error "incorrectly accepted $options"
    }
}

set scratch [file join [pwd] RIVASand_bias_modes_[pid]_[clock clicks]]
file mkdir $scratch
set aliases [list {} {-noBiasVolume} {-fieldBiasVolume}]
foreach {caseName initialStress} {
    dss40 {-19.4 -19.4 -40.0 0.0 0.0 10.0}
    confined60 {-60.0 -60.0 -60.0 0.0 0.0 8.0}
} {
    for {set mode 0} {$mode <= 2} {incr mode} {
        set legacy [runHistory [lindex $aliases $mode] $mode $initialStress ""]
        set selected [runHistory [list -BiasVolume $mode] $mode $initialStress \
            [file join $scratch ${caseName}_${mode}]]
        assertClose $selected $legacy 0.0 "$caseName mode $mode versus previous input"
        set histories($mode) $selected
    }
    if {$caseName eq "confined60"} {
        # Ensure this exercises active differences, not only flags at rest.
        foreach mode {1 2} {
            set maximumStressDifference 0.0
            # Each step stores six stresses followed by 139 state values.
            for {set offset 0} {$offset < [llength $histories(0)]} {incr offset 145} {
                for {set component 0} {$component < 6} {incr component} {
                    set index [expr {$offset+$component}]
                    set difference [expr {abs([lindex $histories($mode) $index]-
                        [lindex $histories(0) $index])}]
                    set maximumStressDifference [expr {max($maximumStressDifference,$difference)}]
                }
            }
            if {$maximumStressDifference <= 1.0e-10} {
                error "mode $mode did not affect the active high-confinement history"
            }
        }
    }
    puts "PASS: $caseName all three modes match previous cyclic histories and restart"
}

# Same-mode duplicates and aliases are accepted in either order.
foreach options {
    {-BiasVolume 0 -BiasVolume 0}
    {-BiasVolume 1 -noBiasVolume}
    {-noBiasVolume -BiasVolume 1}
    {-BiasVolume 2 -fieldBiasVolume}
    {-fieldBiasVolume -BiasVolume 2}
} {
    wipe
    model BasicBuilder -ndm 3 -ndf 3
    {*}[materialCommand $options]
}
wipe
file delete -force $scratch
puts "PASS: BiasVolume parsing, alias equivalence, material copy, and cyclic restart"
