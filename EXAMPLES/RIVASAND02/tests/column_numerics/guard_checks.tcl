# Verify trial-path independence and explicit rejection of unsupported options.
load [lindex $argv 0]
set output [lindex $argv 1]
file mkdir $output
set material RIVASAND02BranchReversalResearch
set row {0.90 900 1.25 1.125 122.44207260468994 0.945 0.025 0.78 0.51 4 4 0.65}
proc makeSpecimen {material row} {
    wipe
    model BasicBuilder -ndm 3 -ndf 4
    nDMaterial $material 1 {*}$row -stage 1 -initialStress -100 -100 -100 0 0 0 -nSub 4 -noBiasVolume
    foreach {n x y z} {1 0 0 0 2 1 0 0 3 1 1 0 4 0 1 0 5 0 0 1 6 1 0 1 7 1 1 1 8 0 1 1} {
        node $n $x $y $z
    }
    element SSPbrickUP 1 1 2 3 4 5 6 7 8 1 550000 1 1e-5 1e-5 1e-5 0.537 2.381e-6 0 0 0
    if {[columnProbe materialIdentity 1] ne $material} {
        error "Material copy lost its research identity"
    }
    columnProbe update
    columnProbe commit
}
proc shearTrial {gamma} {
    foreach n {5 6 7 8} { setNodeDisp $n 1 $gamma }
    columnProbe update
}
model BasicBuilder -ndm 3 -ndf 4
if {![catch {nDMaterial $material 1 {*}$row -reversalLatch}]} {
    error "Research material incorrectly accepted -reversalLatch"
}
set results {}
foreach speculative {0 1} {
    makeSpecimen $material $row
    for {set i 1} {$i<=100} {incr i} {
        shearTrial [expr {$i*.0001}]
        columnProbe commit
    }
    if {$speculative} { shearTrial 0.0101 }
    shearTrial 0.0099
    columnProbe commit
    lappend results [concat [eleResponse 1 stress] [eleResponse 1 state]]
}
foreach a [lindex $results 0] b [lindex $results 1] {
    if {![expr {abs($a-$b)<=1e-13}]} {error "Speculative trial changed accepted state"}
}
database File "$output/disabled_checkpoint"
if {![catch {save 1}]} { error "Research checkpoint unexpectedly succeeded" }
puts "PASS: research latch rejection, trial-path independence, checkpoint rejection"
wipe
