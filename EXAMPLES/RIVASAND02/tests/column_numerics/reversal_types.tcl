# Public selector/copy/trial-state contracts. Requires column_probe.cpp.
load [lindex $argv 0]
set output [lindex $argv 1]
file mkdir $output
set row {0.64 483.48301127222084 1.25 1.125 122.44207260468994 0.945 0.025 0.78 0.51 10 1.5 0.65}

proc specimen {name options expected} {
    global row
    wipe
    model BasicBuilder -ndm 3 -ndf 4
    nDMaterial $name 1 {*}$row -stage 1 -initialStress -100 -100 -100 0 0 0 -nSub 16 -BiasVolume 1 {*}$options
    foreach {n x y z} {1 0 0 0 2 1 0 0 3 1 1 0 4 0 1 0 5 0 0 1 6 1 0 1 7 1 1 1 8 0 1 1} {
        node $n $x $y $z
    }
    element SSPbrickUP 1 1 2 3 4 5 6 7 8 1 550000 1 1e-5 1e-5 1e-5 0.6 2.381e-6 0 0 0
    if {[columnProbe materialIdentity 1] ne $name || [eleResponse 1 reversalType] != $expected} {
        error "Material copy lost identity or reversal type"
    }
    columnProbe update
    columnProbe commit
}

proc setDisplacements {gx gy ez} {
    foreach n {5 6 7 8} {
        setNodeDisp $n 1 $gx -commit
        setNodeDisp $n 2 $gy -commit
        setNodeDisp $n 3 $ez -commit
    }
}
proc trial {gx gy ez} {
    setDisplacements $gx $gy $ez
    columnProbe update
    return [concat [eleResponse 1 stress] [eleResponse 1 state]]
}

proc history {name options expected speculative} {
    specimen $name $options $expected
    set result {}
    set previous {0 0 0}
    for {set i 1} {$i<=256} {incr i} {
        set t [expr {double($i)/128}]
        set gx [expr {.003*sin(2*acos(-1)*$t)}]
        set gy [expr {.0015*(1-cos(2*acos(-1)*$t))}]
        set ez [expr {.0001*$t}]
        set committed [concat [eleResponse 1 stress] [eleResponse 1 state]]
        set reversals [eleResponse 1 reversals]
        if {$speculative} {
            trial [expr {$gx+.002}] [expr {$gy-.001}] [expr {$ez-.0002}]
            # setNodeDisp -commit preserves earlier nodal components but also
            # changes nodal history. Restore it before Domain::revert updates
            # the element again from the reverted material state.
            setDisplacements {*}$previous
            columnProbe revert
            if {[concat [eleResponse 1 stress] [eleResponse 1 state]] ne $committed} {
                error "Type $expected revert did not restore committed state"
            }
            # Also check replacement of an unreverted Newton trial.
            trial [expr {$gx-.002}] [expr {$gy+.001}] [expr {$ez+.0002}]
        }
        lappend result [trial $gx $gy $ez]
        set events [expr {[eleResponse 1 reversals]-$reversals}]
        if {$events<0 || $events>1} {error "More than one reversal in a host increment"}
        columnProbe commit
        set previous [list $gx $gy $ez]
    }
    return $result
}

# Reject invalid, missing, duplicate-conflicting and incompatible selections.
foreach {name options} {
    RIVASAND02 {-reversalType 0}
    RIVASAND02 {-reversalType 4}
    RIVASAND02 {-reversalType -1}
    RIVASAND02 {-reversalType 1.5}
    RIVASAND02 {-reversalType bad}
    RIVASAND02 {-reversalType}
    RIVASAND02 {-reversalType 1 -reversalType 2}
    RIVASAND02 {-reversalType 3 -reversalType 1}
    RIVASAND02 {-reversalType 2 -reversalLatch}
    RIVASAND02 {-reversalLatch -reversalType 2}
    RIVASAND02 {-reversalType 3 -reversalLatch}
    RIVASAND02 {-reversalLatch -reversalType 3}
    RIVASAND02BranchReversalResearch {-reversalType 1}
    RIVASAND02BranchReversalResearch {-reversalType 3}
} {
    wipe
    model BasicBuilder -ndm 3 -ndf 4
    if {![catch {nDMaterial $name 1 {*}$row {*}$options}]} {
        error "Accepted unsupported options: $name $options"
    }
}

foreach mode {1 2 3} {
    set options [list -reversalType $mode -reversalType $mode]
    set baseline($mode) [history RIVASAND02 $options $mode 0]
    foreach tangent {0 1} {
        if {[history RIVASAND02 [concat $options [list -TanType $tangent]] $mode 1] ne $baseline($mode)} {
            error "Type $mode stress/state depends on speculative trials or tangent choice"
        }
    }
    if {$mode!=1} {
        database File "$output/unsupported_checkpoint_$mode"
        if {![catch {save 1}]} {error "Research type $mode checkpoint unexpectedly succeeded"}
    }
}
if {[history RIVASAND02 {} 1 0] ne $baseline(1)} {error "Default type 1 changed"}
if {[history RIVASAND02BranchReversalResearch {} 2 0] ne $baseline(2)} {error "Default alias type 2 changed"}
if {[history RIVASAND02BranchReversalResearch {-reversalType 2} 2 0] ne $baseline(2)} {
    error "Explicit alias type 2 differs"
}
if {$baseline(1) eq $baseline(2) || $baseline(2) eq $baseline(3) || $baseline(1) eq $baseline(3)} {
    error "Rotating path did not distinguish the three rules"
}
specimen RIVASAND02 {-reversalType 1 -reversalLatch} 1
if {[eleResponse 1 reversalLatch]!=1} {error "Type 1 lost latch support"}
puts "REVERSAL_TYPES_COMPLETE: parser, copies, defaults, alias, rollback, trial independence, research guards"
wipe
